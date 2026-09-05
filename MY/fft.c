#include "fft.h"
#include "arm_math.h"
#include <math.h>
#include <string.h>

/*
    主函数流程：读入ADC采样值，存入adc_buf数组，然后对adcbuf数组进行归一化，去直流得到adc_output，
    遍历找最大电压最小电压得到峰峰值，计算有效值.然后对归一化采样数据进行FFT，计算幅值，通过寻峰
    插值等手段找到基频
*/

//==================== 参数配置区 ====================
#define fft_len         4096U        //采样点数
#define Fs              3000000.0f    //采样率 3MHz
#define FREQ_RES 				(Fs/fft_len)
#define FRONT_GAIN      1.0f					//硬件放大倍数
#define VREF_MV         3300.0f				//3.3电压
#define ADC_MAX_CODE    65535.0f				//ADC最高参考值
#define V_FULL  (VREF_MV / 2.0f / FRONT_GAIN) //电压归一化还原系数
#define MAG_THRESHOLD   3.0f         //噪声幅度阈值，低于此值不认为是有效峰
#define PEAK_MIN_MAG    30.0f
#define MAX_PEAK_NUM    60U          //允许峰值数量
#define Freq_Step       500.0f         //频率步进
#define FREQ_TOLERANCE      200.0f      //规整范围
#define CLUSTER_GAP         4000.0f     //聚类阈值
#define FREQ_LOWER      9500.0f     //所有有效分量下限 10kHz
#define FREQ_UPPER      501000.0f    //所有有效分量上限 500kHz
#define HANNING_K           1.5f         //MacLeod汉宁窗修正系数
#define MAX_GROUP_COMPONENT 3U  //最多基波+2谐波
#define FREQ_TOL           3000.0f    //谐波匹配容差
#define MAX_PEROID_NUM      500U
#define chazhi       0.93f
//====================================================

arm_cfft_instance_f32 scfft;          //FFT实例
uint16_t adc_buf[fft_len] = { 0 };					//ADC原始采样数组

float32_t adc_output[fft_len] = { 0 };			//ADC采样数组归一化
float32_t FFT_INPUT[fft_len * 2U] = { 0 }; //FFT复数输入缓存（交错存放Re/Im）
float32_t FFT_OUTPUT[fft_len] = { 0 };
float32_t magSpectrum[fft_len / 2U] = { 0 }; //幅度谱，仅用于寻峰
float32_t window_buffer[fft_len] = { 0 };       //汉宁窗系数缓存
static float32_t window_sum = 0.0f;           //窗函数相干增益的分子 sum(w[n])
float32_t buf_filtered[fft_len];
PeakInfo_t peakBuf[MAX_PEAK_NUM];
uint16_t peakCnt = 0;;                   //合法峰值个数
PeakInfo_t finalPeak[MAX_PEAK_NUM];   //分簇+规整之后输出峰
uint8_t peakOverflowFlag = 0;


// 生成汉宁窗系数，main中调用一次
void Generate_hanning_window() {
    window_sum = 0.0f;
    for (uint16_t i = 0; i < fft_len; i++) {
        window_buffer[i] = 0.5f * (1.0f - cosf(2.0f * PI * i / (fft_len - 1)));
        window_sum += window_buffer[i];
    }
}

//初始化
void FFT_Init() {
    arm_cfft_init_f32(&scfft, fft_len);
    memset(adc_output, 0, sizeof(adc_output));
    memset(FFT_INPUT, 0, sizeof(FFT_INPUT));
    memset(magSpectrum, 0, sizeof(magSpectrum));
    Generate_hanning_window();
}



//ADC数据转换
void ADC_ArrayToVoltage(uint16_t* adc_in, float32_t* adc_out, uint16_t len) {
    float32_t mean = 0.0f;
    for (uint16_t i = 0; i < len; i++) {
        //先归一化
        adc_out[i] = (adc_in[i] / 32767.5f) - 1.0f;
        //adc_out[i] = (adc_in[i] / 65525.0f * 2.0f);
        //递推滑动均值，防止大数累积误差
        float32_t val = adc_out[i];
        mean = mean + (val - mean) / (i + 1.0f);
    }
    //时域每个采样减去直流，得到交流信号
    for (uint16_t i = 0; i < len; i++) {
        adc_out[i] = adc_out[i] - mean;
    }
}


//FFT处理
void FFT_Process(float32_t* fft_in) {
    //时域加汉宁窗，填入FFT实部，虚部置0
    for (uint16_t i = 0; i < fft_len; i++) {
        FFT_INPUT[2 * i] = fft_in[i] * window_buffer[i];
        FFT_INPUT[2 * i + 1] = 0.0f;
    }
    arm_cfft_f32(&scfft, FFT_INPUT, 0, 1);
    arm_cmplx_mag_f32(FFT_INPUT, magSpectrum, fft_len / 2U);

}

/*
//插值
static uint8_t MacLeodInterpLegacy(float32_t *fftComplex, uint16_t k0, PeakInfo_t *peakOut)
{
    const uint16_t halfN = fft_len / 2U;
    //最左/最右频点缺少一侧邻点，无法三点插值，降级处理
    if((k0 == 0U) || (k0 >= halfN - 1U))
    {
        peakOut->freq = (float32_t)k0 * Fs / fft_len;
        float32_t re = fftComplex[2U * k0];
        float32_t im = fftComplex[2U * k0 + 1U];
        peakOut->mag = sqrtf(re * re + im * im);
        peakOut->phase = atan2f(im, re);
        return 1U;
    }

    //读取k0-1 左邻点复数
    float32_t Xm_re = fftComplex[2U * (k0 - 1U)];
    float32_t Xm_im = fftComplex[2U * (k0 - 1U) + 1U];
    //读取当前k0复数
    float32_t X0_re = fftComplex[2U * k0];
    float32_t X0_im = fftComplex[2U * k0 + 1U];
    //读取k0+1 右邻点复数
    float32_t Xp_re = fftComplex[2U * (k0 + 1U)];
    float32_t Xp_im = fftComplex[2U * (k0 + 1U) + 1U];

    float32_t num_re = Xp_re - Xm_re;
    float32_t num_im = Xp_im - Xm_im;
    float32_t den_re = 2.0f * X0_re - Xp_re - Xm_re;
    float32_t den_im = 2.0f * X0_im - Xp_im - Xm_im;

    float32_t den_sq = den_re * den_re + den_im * den_im;
    //分母趋近0，栅格正中，无偏移，降级
    if(fabsf(den_sq) < 1e-12f)
    {
        peakOut->freq = (float32_t)k0 * Fs / fft_len;
        peakOut->mag = sqrtf(X0_re * X0_re + X0_im * X0_im);
        peakOut->phase = atan2f(X0_im, X0_re);
        return 1U;
    }

    //求解栅格偏移量 delta
    float32_t div_re = (num_re * den_re + num_im * den_im) / den_sq;
    float32_t delta_raw = 0.5f * div_re;
    float32_t delta = HANNING_K * delta_raw;

    //浮点等效栅格索引，换算真实频率
    float32_t k_hat = (float32_t)k0 + delta;
    peakOut->freq = k_hat * Fs / fft_len;

    //汉宁窗频域幅度补偿，还原真实信号幅值
    float32_t pi_d = PI * delta;
    float32_t Wd;
    if(fabsf(delta) < 1e-10f)
    {
        Wd = 1.0f;
    }
    else
    {
        Wd = sinf(pi_d) / (pi_d * (1.0f - delta * delta));
    }
    Wd = fabsf(Wd);

    float32_t mag0 = sqrtf(X0_re * X0_re + X0_im * X0_im);
    peakOut->mag = mag0 / Wd;
    peakOut->phase = atan2f(X0_im, X0_re);

    return 0U;
}
*/


#if 0
/* 旧实现仅留作对照：混用了三种估计器及经验幅值系数。 */
static uint8_t MacLeodInterp(float32_t* fftComplex, uint16_t k0, PeakInfo_t* peakOut) {
    const uint16_t halfN = fft_len / 2U;

    // 最左/最右频点缺少一侧邻点，无法三点插值，降级处理
    if ((k0 == 0U) || (k0 >= halfN - 1U)) {
        peakOut->freq = (float32_t)k0 * Fs / fft_len;
        float32_t re = fftComplex[2U * k0];
        float32_t im = fftComplex[2U * k0 + 1U];
        peakOut->mag = sqrtf(re * re + im * im);
        peakOut->phase = atan2f(im, re);
        return 1U;
    }

    // 读取三点复数
    float32_t Xm_re = fftComplex[2U * (k0 - 1U)];
    float32_t Xm_im = fftComplex[2U * (k0 - 1U) + 1U];
    float32_t X0_re = fftComplex[2U * k0];
    float32_t X0_im = fftComplex[2U * k0 + 1U];
    float32_t Xp_re = fftComplex[2U * (k0 + 1U)];
    float32_t Xp_im = fftComplex[2U * (k0 + 1U) + 1U];

    // 计算三点幅度
    float32_t mag_m = sqrtf(Xm_re * Xm_re + Xm_im * Xm_im);
    float32_t mag_0 = sqrtf(X0_re * X0_re + X0_im * X0_im);
    float32_t mag_p = sqrtf(Xp_re * Xp_re + Xp_im * Xp_im);

    // 防止除零
    if (mag_0 < 1e-12f) {
        peakOut->freq = (float32_t)k0 * Fs / fft_len;
        peakOut->mag = 0.0f;
        peakOut->phase = 0.0f;
        return 1U;
    }

    // ============ 第一步：频率偏移量 delta 计算 ============
    float32_t delta = 0;

    // 使用幅度比法（对频谱泄露更鲁棒）
    float32_t r1 = mag_m / mag_0;
    float32_t r2 = mag_p / mag_0;

    float32_t delta_num = (r1 - r2) * 2.0f;
    float32_t delta_den = (r1 + r2) * 2.0f - 4.0f;

    if (fabsf(delta_den) < 1e-12f) {
        delta = 0.0f;
    }
    else {
        delta = delta_num / delta_den;
        // 汉宁窗校正
        delta = delta * HANNING_K;
    }

    // ============ 第二步：二次精细化（提升半格精度） ============
    // 当偏移量接近 ±0.5 时，使用二次插值修正
    if (fabsf(delta) > 0.35f) {
        // 方法2：抛物线插值（作为参考）
        float32_t delta_para = 0.5f * (mag_m - mag_p) / (mag_m - 2.0f * mag_0 + mag_p);
        delta_para = delta_para * HANNING_K;

        // 方法3：复数比值法（作为参考）
        float32_t num_re = Xp_re - Xm_re;
        float32_t num_im = Xp_im - Xm_im;
        float32_t den_re = 2.0f * X0_re - Xp_re - Xm_re;
        float32_t den_im = 2.0f * X0_im - Xp_im - Xm_im;
        float32_t den_sq = den_re * den_re + den_im * den_im;
        float32_t delta_cplx = 0.0f;
        if (fabsf(den_sq) > 1e-12f) {
            float32_t div_re = (num_re * den_re + num_im * den_im) / den_sq;
            delta_cplx = 0.5f * div_re * HANNING_K;
        }

        // 加权平均：幅度比法权重0.5，抛物线0.25，复数法0.25
        delta = 0.5f * delta + 0.25f * delta_para + 0.25f * delta_cplx;
    }

    // 限幅，防止发散
    if (delta > 0.5f) delta = 0.5f;
    if (delta < -0.5f) delta = -0.5f;

    // ============ 第三步：频率计算 ============
    float32_t k_hat = (float32_t)k0 + delta;
    peakOut->freq = k_hat * Fs / fft_len;

    // ============ 第四步：幅度恢复 ============
    // 确定积分范围
    //汉宁窗频域幅度补偿，还原真实信号幅值
    float32_t pi_d = PI * delta;
    float32_t Wd = 0;
    if (fabsf(delta) < 1e-10f) {
        Wd = 1.0f;
    }
    else {
        Wd = sinf(pi_d) / (pi_d * (1.0f - delta * delta));
    }
    Wd = fabsf(Wd);

    float32_t mag0 = sqrtf(X0_re * X0_re + X0_im * X0_im);
    peakOut->mag = mag0 / Wd * 0.91f;
    // ============ 第五步：相位计算（带校正） ============
    // 基础相位
    float32_t phase0 = atan2f(X0_im, X0_re);

    // 相位校正（消除偏移带来的相位误差）
    float32_t phase_correction = PI * delta * (fft_len - 1.0f) / fft_len;
    peakOut->phase = phase0 + phase_correction;

    // 归一化到 [-PI, PI]
    while (peakOut->phase > PI) peakOut->phase -= 2.0f * PI;
    while (peakOut->phase < -PI) peakOut->phase += 2.0f * PI;

    return 0U;
}
#endif

/*
 * 汉宁窗三点频域插值。
 * 输出 mag 已经是输入时域正弦量的峰值，不再保留原始 FFT 标度。
 */
static uint8_t HannPeakInterp(const float32_t* fftComplex,
    uint16_t k0,
    PeakInfo_t* peakOut) {
    const uint16_t halfN = fft_len / 2U;
    const float32_t eps = 1.0e-12f;
    float32_t re0 = fftComplex[2U * k0];
    float32_t im0 = fftComplex[2U * k0 + 1U];
    float32_t mag0 = sqrtf(re0 * re0 + im0 * im0);
    float32_t delta = 0.0f;
    uint8_t degraded = 0U;

    if ((k0 == 0U) || (k0 >= halfN - 1U) || (mag0 <= eps)) {
        degraded = 1U;
    }
    else {
        float32_t rem = fftComplex[2U * (k0 - 1U)];
        float32_t imm = fftComplex[2U * (k0 - 1U) + 1U];
        float32_t rep = fftComplex[2U * (k0 + 1U)];
        float32_t imp = fftComplex[2U * (k0 + 1U) + 1U];
        float32_t magm = sqrtf(rem * rem + imm * imm);
        float32_t magp = sqrtf(rep * rep + imp * imp);
        float32_t den = magm + 2.0f * mag0 + magp;

        /* 汉宁窗专用三点估计式：正值表示峰值位于 k0 右侧。 */
        if (den > eps)
            delta = 2.0f * (magp - magm) / den;

        if (delta > 0.5f)  delta = 0.5f;
        if (delta < -0.5f) delta = -0.5f;
    }

    peakOut->freq = ((float32_t)k0 + delta) * FREQ_RES;

    /*
     * 汉宁窗主瓣在 delta 处的归一化响应：
     * sinc(delta)/(1-delta^2)。
     * 实信号单边谱峰值 A = 2*|X[k0]|/(sum(window)*response)。
     */
    {
        float32_t response = 1.0f;
        if (fabsf(delta) > 1.0e-6f) {
            float32_t x = PI * delta;
            response = fabsf(sinf(x) / (x * (1.0f - delta * delta)));
        }
        peakOut->mag = (window_sum > eps && response > eps)
            ? (2.0f * mag0 / (window_sum * response) * chazhi)
            : 0.0f;
    }

    peakOut->phase = atan2f(im0, re0)
        + PI * delta * (fft_len - 1.0f) / fft_len;
    while (peakOut->phase > PI)  peakOut->phase -= 2.0f * PI;
    while (peakOut->phase < -PI) peakOut->phase += 2.0f * PI;

    return degraded;
}


static float32_t RoundFreqTo500Grid(float32_t freq) {
    float32_t kf = freq / 500.0f;
    int32_t k_round = (int32_t)(kf + 0.5f);
    float32_t f_round = (float32_t)k_round * 500.0f;
    if (fabsf(freq - f_round) < 250.0f) {
        return f_round;
    }
    return freq;
}

// 峰值数组按照频率从小到大冒泡排序
void SortPeakByFreq(PeakInfo_t* p, uint16_t num) {
    for (uint16_t i = 0; i < num; i++) {
        for (uint16_t j = i + 1; j < num; j++) {
            if (p[j].freq < p[i].freq) {
                PeakInfo_t tmp = p[i];
                p[i] = p[j];
                p[j] = tmp;
            }
        }
    }
}


//泄露簇聚类清洗
static uint16_t RemoveLeakCluster(PeakInfo_t* inPeak, uint16_t inNum, PeakInfo_t* outPeak) {
    if (inNum == 0)
        return 0;

    SortPeakByFreq(inPeak, inNum);
    uint16_t outCnt = 0;
    uint16_t startIdx = 0;

    while (startIdx < inNum) {
        uint16_t endIdx = startIdx;
        //向后扩张，寻找整个连续簇
        while ((endIdx + 1 < inNum) &&
            fabsf(inPeak[endIdx + 1].freq - inPeak[endIdx].freq) < CLUSTER_GAP) {
            endIdx++;
        }
        //簇内寻找幅值最大点作为真实主峰
        uint16_t maxIdx = startIdx;
        for (uint16_t i = startIdx; i <= endIdx; i++) {
            if (inPeak[i].mag > inPeak[maxIdx].mag)
                maxIdx = i;
        }
        outPeak[outCnt++] = inPeak[maxIdx];
        startIdx = endIdx + 1U;
    }
    return outCnt;
}

//寻峰
static void SpectrumFindPeak(float32_t* mag, float32_t* fftComplex) {
    PeakInfo_t tempPeakBuf[MAX_PEAK_NUM];
    uint16_t tempPeakCnt = 0;
    PeakInfo_t clusterBuf[MAX_PEAK_NUM];
    uint16_t clusterCnt = 0;
    peakCnt = 0;
    peakOverflowFlag = 0; //每次寻峰清空溢出标记

    const uint16_t halfN = fft_len / 2U;
    uint16_t k_max = (uint16_t)(FREQ_UPPER / (Fs / fft_len));
    if (k_max >= halfN)
        k_max = halfN - 1U;

    for (uint16_t k = 1; k < k_max; k++) {
        //仅判断左右相邻1点，噪声下依然可以识别峰，放宽要求避免平顶峰丢失主峰
        if ((mag[k] >= mag[k - 1]) && (mag[k] >= mag[k + 1]) &&
            !((mag[k] == mag[k - 1]) && (mag[k] == mag[k + 1]))) {
            if (mag[k] > MAG_THRESHOLD) {
                if (tempPeakCnt < MAX_PEAK_NUM) {
                    PeakInfo_t tmpPeak;
                    HannPeakInterp(fftComplex, k, &tmpPeak);
                    tempPeakBuf[tempPeakCnt++] = tmpPeak;
                }
                else {
                    peakOverflowFlag = 1; //缓存已满，部分峰值丢弃
                }
            }
        }
    }

    clusterCnt = RemoveLeakCluster(tempPeakBuf, tempPeakCnt, clusterBuf);
    // 直接输出到peakBuf
    peakCnt = clusterCnt;
    memcpy(peakBuf, clusterBuf, peakCnt * sizeof(PeakInfo_t));
}

//谐波匹配
HarmonicGroup_t signalGroup;
void HarmonicMatch(PeakInfo_t* peakBuf, uint16_t peakCnt, HarmonicGroup_t* bestGroup) {
    //初始化最优组
    bestGroup->cnt = 0;
    bestGroup->sumMag = 0.0f;

    //全部峰值标记为未占用
    for (uint16_t i = 0; i < peakCnt; i++) {
        peakBuf[i].used = 0;
    }

    //遍历每个峰值作为候选基波
    for (uint16_t i = 0; i < peakCnt; i++) {
        if (peakBuf[i].used != 0)
            continue;

        float32_t f0 = peakBuf[i].freq;
        //频域边界过滤
        if (f0 < FREQ_LOWER || f0 > FREQ_UPPER)
            continue;

        HarmonicGroup_t tempGroup;
        tempGroup.cnt = 0;
        tempGroup.sumMag = 0.0f;

        //基波加入临时组
        tempGroup.peaks[tempGroup.cnt++] = peakBuf[i];
        tempGroup.sumMag += peakBuf[i].mag;

        //向后搜寻谐波 j>i，天然保证基波是组内最低频率，杜绝反向匹配
        for (uint16_t j = i + 1; j < peakCnt; j++) {
            if (peakBuf[j].used != 0)
                continue;

            float32_t fh = peakBuf[j].freq;
            if (fh < FREQ_LOWER || fh > FREQ_UPPER)
                continue;

            //计算理论谐波次数，就近取整
            float32_t n_float = fh / f0;
            uint16_t n = (uint16_t)roundf(n_float);
            if (n < 2)
                continue;   //谐波次数≥2

            float32_t f_theory = (float32_t)n * f0;
            //频率在容差范围内，判定匹配成功
            if (fabsf(fh - f_theory) < FREQ_TOL) {
                tempGroup.peaks[tempGroup.cnt++] = peakBuf[j];
                tempGroup.sumMag += peakBuf[j].mag;
                if (tempGroup.cnt >= MAX_GROUP_COMPONENT)
                    break; //达到最大数量，停止搜索更高次谐波
            }
        }

        //择优更新最优组
        uint8_t update_flag = 0;
        if (tempGroup.cnt > bestGroup->cnt) {
            update_flag = 1;
        }
        else if (tempGroup.cnt == bestGroup->cnt && tempGroup.sumMag > bestGroup->sumMag) {
            update_flag = 1;
        }

        if (update_flag) {
            *bestGroup = tempGroup;
        }
    }
}



//主函数
uint8_t Signal_Process(float32_t* Vpp, float32_t* Vrms, float32_t* fA, float32_t* FREQ, float32_t* Mag) {
    (void)Vpp;
    *Vrms = 0.0f;
    *fA = 0.0f;
    for (uint8_t i = 0; i < MAX_GROUP_COMPONENT; i++) {
        FREQ[i] = 0.0f;
        Mag[i] = 0.0f;
    }

    //初始化
    FFT_Init();

    //ADC数据转换
    ADC_ArrayToVoltage(adc_buf, adc_output, fft_len);

    //进行FFT计算
    FFT_Process(adc_output);

    //寻峰
    SpectrumFindPeak(magSpectrum, FFT_INPUT);
    SortPeakByFreq(peakBuf, peakCnt);

    if (peakCnt == 0U) {
        signalGroup.cnt = 0U;
        signalGroup.sumMag = 0.0f;
        return 0U;
    }

    *fA = RoundFreqTo500Grid(peakBuf[0].freq);
    //谐波匹配
    HarmonicMatch(peakBuf, peakCnt, &signalGroup);
    if (signalGroup.cnt >= 1) {
        for (uint8_t i = 0; i < signalGroup.cnt; i++) {
            FREQ[i] = RoundFreqTo500Grid(signalGroup.peaks[i].freq);
            /* 插值函数已经完成单边谱及汉宁窗增益归一化。 */
            Mag[i] = signalGroup.peaks[i].mag * V_FULL;
            *Vrms += Mag[i] * Mag[i];
        }
        *Vrms = sqrtf(*Vrms / 2.0f);

    }
    for (uint8_t i = 1; i < signalGroup.cnt; i++) {
        // 计算理论谐波次数，四舍五入取整
        float ratio = FREQ[i] / FREQ[0];
        uint16_t n = (uint16_t)(ratio + 0.5f);  // 四舍五入到整数

        // 防错：次数不能为0，至少1次（基波自身）
        if (n < 1) n = 1;

        // 规整后的频率 = 基波 × 整数次数
        FREQ[i] = FREQ[0] * n;
    }
    return (signalGroup.cnt > 0U) ? 1U : 0U;
}
