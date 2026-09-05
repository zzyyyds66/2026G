#include "fft8192.h"
#include <math.h>
#include <string.h>

#define PI  3.14159265358979f

// 缓冲区定义
fft32_t fft8192_real[FFT8192_LEN] = { 0 };          //FFT实部缓冲
fft32_t fft8192_imag[FFT8192_LEN] = { 0 };          //FFT虚部缓冲
fft32_t fft8192_mag[FFT8192_LEN / 2U] = { 0 };       //幅度谱(仅前N/2点有效)
fft32_t fft8192_window[FFT8192_LEN] = { 0 };         //汉宁窗系数
fft32_t fft8192_twiddle_cos[FFT8192_LEN / 2U] = { 0 }; //旋转因子cos表
fft32_t fft8192_twiddle_sin[FFT8192_LEN / 2U] = { 0 }; //旋转因子sin表
fft32_t fft8192_phase[FFT8192_LEN / 2U] = { 0 };      //相位谱
FFT8192_Peak_t fft8192_peaks[FFT8192_MAX_PEAKS];      //峰值数组
uint16_t fft8192_peakCnt = 0;                          //峰值数量
fft32_t fft8192_magcorr[FFT8192_LEN / 2U] = { 0 };       //插值与补偿之后的幅度谱(仅前N/2点有效)

fft32_t window_sum = 0.0f;                    //窗函数相干增益
fft32_t fft8192_vrms = 0.0f;                  //真有效值(mV)
fft32_t fft8192_vpp = 0.0f;                   //峰峰值(mV)
fft32_t fft8192_out_freq[3] = { 0 };          //输出频率数组(Hz)
fft32_t fft8192_out_mag[3] = { 0 };           //输出幅度数组(mV)
fft32_t fft8192_out_phase[3] = { 0 };         //输出相位数组(rad)


// 输入滤波器增益补偿表(10kHz~500kHz, 实测1/gain)
#define FILTER_TABLE_LEN  50
fft32_t filter_gain_table[FILTER_TABLE_LEN] = {
    0.9363f, 0.9355f, 0.9420f, 0.9440f, 0.9440f,
    0.9440f, 0.9440f, 0.9440f, 0.9440f, 0.9420f,
    0.9390f, 0.9410f, 0.9340f, 0.9380f, 0.9460f,
    0.9380f, 0.9390f, 0.9400f, 0.9410f, 0.9430f,
    0.9420f, 0.9440f, 0.9470f, 0.9500f, 0.9510f,
    0.9550f, 0.9655f, 0.9670f, 0.9680f, 0.9710f,
    0.9710f, 0.9725f, 0.9740f, 0.9760f, 0.9810f,
    0.9850f, 0.9880f, 0.9890f, 0.9920f, 0.9980f,
    1.0025f, 1.0085f, 1.0120f, 1.0230f, 1.0320f,
    1.0400f, 1.0500f, 1.0630f, 1.0750f, 1.0900f,
};
#define FILTER_FREQ_START  10000.0f
#define FILTER_FREQ_STEP   10000.0f

// 输入滤波器相位延迟表(10kHz~500kHz, 实测度数, 已解卷绕)
fft32_t filter_phase_table[FILTER_TABLE_LEN] = {
    5.2f, 9.0f, 11.9f, 14.5f, 18.5f,
    21.5f, 27.5f, 30.25f, 35.25f, 38.5f,
    43.5f, 47.5f, 49.25f, 53.0f, 57.25f,
    63.0f, 67.25f, 68.0f, 72.0f, 76.5f,
    82.0f, 86.0f, 88.75f, 94.75f, 98.75f,
    102.0f, 105.5f, 111.0f, 115.75f, 120.0f,
    123.5f, 127.5f, 132.0f, 136.5f, 140.5f,
    146.0f, 150.5f, 154.0f, 157.5f, 163.5f,
    167.5f, 174.0f, 177.5f, 183.5f, 188.5f,
    193.5f, 201.5f, 205.5f, 210.5f, 217.0f,
};

// 位反转排列
static void BitReverse(fft32_t* re, fft32_t* im) {
    uint16_t i, j = 0, k;
    for (i = 0; i < FFT8192_LEN - 1; i++) {
        if (i < j) {
            fft32_t tr = re[i]; re[i] = re[j]; re[j] = tr;
            fft32_t ti = im[i]; im[i] = im[j]; im[j] = ti;
        }
        k = FFT8192_LEN >> 1;
        while (k <= j) { j -= k; k >>= 1; }
        j += k;
    }
}

// 预计算旋转因子 W_N^k = cos(-2πk/N) + j*sin(-2πk/N)
static void GenerateTwiddle(void) {
    for (uint16_t k = 0; k < FFT8192_LEN / 2U; k++) {
        fft8192_twiddle_cos[k] = cosf(2.0f * PI * k / FFT8192_LEN);
        fft8192_twiddle_sin[k] = -sinf(2.0f * PI * k / FFT8192_LEN);
    }
}

// 基2蝶形FFT
static void FFT_Radix2(fft32_t* re, fft32_t* im) {
    for (uint16_t s = 1; s <= FFT8192_LOG2N; s++) {
        uint16_t m = 1U << s;              //当前DFT大小
        uint16_t m2 = m >> 1;              //蝶形两半间距
        uint16_t step = FFT8192_LEN / m;   //旋转因子步进

        for (uint16_t k = 0; k < FFT8192_LEN; k += m) {
            for (uint16_t j = 0; j < m2; j++) {
                uint16_t tw_idx = j * step;
                fft32_t wr = fft8192_twiddle_cos[tw_idx];
                fft32_t wi = fft8192_twiddle_sin[tw_idx];

                uint16_t i = k + j;
                uint16_t ip = i + m2;

                fft32_t tr = re[ip] * wr - im[ip] * wi;
                fft32_t ti = re[ip] * wi + im[ip] * wr;

                re[ip] = re[i] - tr;
                im[ip] = im[i] - ti;
                re[i] += tr;
                im[i] += ti;
            }
        }
    }
}

// 生成汉宁窗系数
void FFT8192_GenerateWindow(void) {
    window_sum = 0.0f;
    for (uint16_t i = 0; i < FFT8192_LEN; i++) {
        fft8192_window[i] = 0.5f * (1.0f - cosf(2.0f * PI * i / (FFT8192_LEN - 1)));
        window_sum += fft8192_window[i];
    }
}

// 初始化: 预计算旋转因子和汉宁窗
void FFT8192_Init(void) {
    memset(fft8192_real, 0, sizeof(fft8192_real));
    memset(fft8192_imag, 0, sizeof(fft8192_imag));
    memset(fft8192_mag, 0, sizeof(fft8192_mag));
    GenerateTwiddle();
    FFT8192_GenerateWindow();
}

// FFT处理: 归一化去直流→加窗→位反转→蝶形运算→求幅度谱和相位
void FFT8192_Process(const uint16_t* adc_in) {
    // 1. 归一化到[-1,1]并递推求均值
    fft32_t mean = 0.0f;
    for (uint16_t i = 0; i < FFT8192_LEN; i++) {
        fft8192_real[i] = ((fft32_t)adc_in[i] / 32767.5f) - 1.0f;
        mean += (fft8192_real[i] - mean) / (fft32_t)(i + 1);
    }
    // 2. 去直流, 计算时域Vrms/Vpp, 加汉宁窗, 虚部置0
    fft32_t sig_max = -1e30f, sig_min = 1e30f, sum_sq = 0.0f;
    for (uint16_t i = 0; i < FFT8192_LEN; i++) {
        fft8192_real[i] -= mean;   // 去直流
        if (fft8192_real[i] > sig_max) sig_max = fft8192_real[i];
        if (fft8192_real[i] < sig_min) sig_min = fft8192_real[i];
        sum_sq += fft8192_real[i] * fft8192_real[i];
        fft8192_real[i] *= fft8192_window[i];   // 加窗
        fft8192_imag[i] = 0.0f;
    }
    fft8192_vpp = (sig_max - sig_min) * 3300.0f;   // 峰峰值(mV)
    fft8192_vrms = sqrtf(sum_sq / FFT8192_LEN) * 3300.0f;   // 真有效值(mV)

    BitReverse(fft8192_real, fft8192_imag);
    FFT_Radix2(fft8192_real, fft8192_imag);

    // 1. 求原始幅度谱和相位(不做校正, 保留原始值供插值使用)
    for (uint16_t i = 0; i < FFT8192_LEN / 2U; i++) {
        fft32_t re = fft8192_real[i];
        fft32_t im = fft8192_imag[i];
        fft32_t val = re * re + im * im;
        if (val < 0.0f) val = 0.0f;
        fft8192_mag[i] = sqrtf(val);
        fft8192_phase[i] = atan2f(im, re);
    }

    // 2. 三谱线插值寻峰(算mag_corr时顺便存入fft8192_magcorr)
    FFT8192_FindPeaks();

    // 3. 谐波匹配: 找基波+谐波, 频率规整为500Hz整数倍
    FFT8192_HarmonicMatch();

    // 4. 复制结果到全局输出数组供其他模块使用
    for (uint8_t i = 0; i < 3; i++) {
        if (i < fft8192_peakCnt) {
            fft8192_out_freq[i] = fft8192_peaks[i].freq;
            fft8192_out_mag[i] = fft8192_peaks[i].mag;
            fft8192_out_phase[i] = fft8192_peaks[i].phase;
        }
        else {
            fft8192_out_freq[i] = 0.0f;
            fft8192_out_mag[i] = 0.0f;
            fft8192_out_phase[i] = 0.0f;
        }
    }
}


// 滤波器增益补偿: 10kHz以下返回1.0, 500kHz以上用末端值, 中间线性插值
fft32_t FFT8192_FilterGain(fft32_t freq) {
    if (freq < FILTER_FREQ_START) return filter_gain_table[0];
    if (freq > 500000.0f) return filter_gain_table[FILTER_TABLE_LEN - 1];
    float idx_f = (freq - FILTER_FREQ_START) / FILTER_FREQ_STEP;
    uint16_t idx = (uint16_t)idx_f;
    float ratio = idx_f - idx;
    return filter_gain_table[idx] + ratio * (filter_gain_table[idx + 1] - filter_gain_table[idx]);
}

// 滤波器相位延迟(弧度): 10kHz以下用首值, 500kHz以上用末值, 中间线性插值
fft32_t FFT8192_PhaseDelay(fft32_t freq) {
    if (freq < FILTER_FREQ_START) return filter_phase_table[0] * PI / 180.0f;
    if (freq > 500000.0f) return filter_phase_table[FILTER_TABLE_LEN - 1] * PI / 180.0f;
    float idx_f = (freq - FILTER_FREQ_START) / FILTER_FREQ_STEP;
    uint16_t idx = (uint16_t)idx_f;
    float ratio = idx_f - idx;
    float deg = filter_phase_table[idx] + ratio * (filter_phase_table[idx + 1] - filter_phase_table[idx]);
    return deg * PI / 180.0f;
}

// 三谱线插值寻峰: 用k-1,k,k+1三根谱线精估计频率/幅度/相位
void FFT8192_FindPeaks(void) {
    memset(fft8192_peaks, 0, sizeof(fft8192_peaks));   //清空峰值数组
    fft8192_peakCnt = 0;
    uint16_t k_max = (uint16_t)(FFT8192_FREQ_UPPER / FFT8192_FREQ_RES);
    if (k_max >= FFT8192_LEN / 2U) k_max = FFT8192_LEN / 2U - 1U;

    for (uint16_t k = 1; k < k_max; k++) {
        if (fft8192_mag[k] < FFT8192_MAG_TH) continue;                                       // 噪声阈值过滤
        if (!(fft8192_mag[k] >= fft8192_mag[k - 1] && fft8192_mag[k] >= fft8192_mag[k + 1])) continue; // 局部极大值判定
        if (fft8192_peakCnt >= FFT8192_MAX_PEAKS) break;                                     // 峰值缓冲上限
        if ((fft32_t)k * FFT8192_FREQ_RES < FFT8192_FREQ_LOWER) continue;                    // 频率下限过滤

        fft32_t mag_m = fft8192_mag[k - 1];   // 左邻谱线幅度
        fft32_t mag_0 = fft8192_mag[k];       // 中心谱线幅度
        fft32_t mag_p = fft8192_mag[k + 1];   // 右邻谱线幅度

        fft32_t den = mag_m + 2.0f * mag_0 + mag_p;   // 三谱线加权和
        fft32_t delta = 0.0f;                          // 频偏量(-0.5~0.5)
        if (den > 1e-12f) delta = 2.0f * (mag_p - mag_m) / den;
        if (delta > 0.5f) delta = 0.5f;
        if (delta < -0.5f) delta = -0.5f;

        fft32_t freq = ((fft32_t)k + delta) * FFT8192_FREQ_RES;   // 精确频率

        fft32_t response = 1.0f;   // 汉宁窗主瓣响应
        if (fabsf(delta) > 1e-6f) {
            fft32_t x = PI * delta;
            response = fabsf(sinf(x) / (x * (1.0f - delta * delta)));
        }
        fft32_t mag_corr = (response > 1e-6f) ? (2.0f * mag_0 / (window_sum * response)) : (2.0f * mag_0 / window_sum);   // 单边谱×2+窗增益+主瓣响应校正
        mag_corr *= FFT8192_FilterGain(freq);   // 滤波器增益补偿
        mag_corr *= 3300.0f;    // ADC量程校正(mV)
        fft8192_magcorr[k] = mag_corr;   // 存储补偿后幅度供其他模块使用
        if (mag_corr < FFT8192_PEAK_MIN) continue;   // 幅度低于阈值视为噪声

        fft32_t phase = fft8192_phase[k] + PI * delta * (fft32_t)(FFT8192_LEN - 1) / FFT8192_LEN;   // 相位校正
        while (phase > PI) phase -= 2.0f * PI;
        while (phase < -PI) phase += 2.0f * PI;

        fft8192_peaks[fft8192_peakCnt].freq = freq;
        fft8192_peaks[fft8192_peakCnt].mag = mag_corr;
        fft8192_peaks[fft8192_peakCnt].phase = phase;
        fft8192_peakCnt++;
    }
}

// 谐波匹配: 基波为频率最低的峰, 搜索其整数倍谐波, 频率规整为500Hz整数倍
void FFT8192_HarmonicMatch(void) {
    if (fft8192_peakCnt == 0) return;

    // 1. 频率最低的峰作为基波, 规整到500Hz栅格
    fft32_t f0_raw = fft8192_peaks[0].freq;
    int32_t grid_n = (int32_t)(f0_raw / FFT8192_FREQ_GRID + 0.5f);
    fft32_t f0 = (fft32_t)grid_n * FFT8192_FREQ_GRID;
    if (f0 < FFT8192_FREQ_LOWER || f0 > FFT8192_FREQ_UPPER) return;

    // 3. 搜索谐波: 频率约为基波整数倍
    FFT8192_Peak_t result[3];
    uint8_t result_cnt = 0;

    result[0].freq = f0;
    result[0].mag = fft8192_peaks[0].mag;
    result[0].phase = fft8192_peaks[0].phase;
    result_cnt = 1;

    for (uint16_t i = 1; i < fft8192_peakCnt && result_cnt < 3; i++) // 搜索谐波，条件为
    {
        float ratio = fft8192_peaks[i].freq / f0; // 谐波次数
        uint16_t n = (uint16_t)(ratio + 0.5f);   // 四舍五入取谐波次数
        if (n < 2) continue;                      // 谐波次数≥2

        fft32_t f_theory = (fft32_t)n * f0;       // 理论谐波频率
        if (fabsf(fft8192_peaks[i].freq - f_theory) > FFT8192_FREQ_GRID) continue;   // 容差±500Hz

        result[result_cnt].freq = f_theory;       // 规整为基波整数倍
        result[result_cnt].mag = fft8192_peaks[i].mag;
        result[result_cnt].phase = fft8192_peaks[i].phase;
        result_cnt++;
    }

    // 4. 回写, 只保留基波+匹配的谐波
    for (uint8_t i = 0; i < result_cnt; i++) {
        fft8192_peaks[i] = result[i];
    }
    fft8192_peakCnt = result_cnt;
}
