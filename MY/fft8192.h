#ifndef FFT8192_H
#define FFT8192_H

#include <stdint.h>

#define FFT8192_LEN     8192U           //采样点数(2^13)
#define FFT8192_FS      3000000.0f      //采样率 3MHz
#define FFT8192_FREQ_RES    (FFT8192_FS / FFT8192_LEN)
#define FFT8192_LOG2N   13U             //log2(8192)

typedef float fft32_t;

// 峰值信息
typedef struct {
    fft32_t freq;       //频率
    fft32_t mag;         //幅度
    fft32_t phase;       //相位
} FFT8192_Peak_t;

#define FFT8192_MAX_PEAKS   5U         //最大基波/谐波峰值数
#define FFT8192_MAG_TH      0.001f      //噪声幅度阈值
#define FFT8192_PEAK_MIN    4.0f       //峰值有效幅度下限(mV)
#define FFT8192_FREQ_LOWER  9500.0f    //有效频率下限 10kHz
#define FFT8192_FREQ_UPPER  501000.0f   //有效频率上限 500kHz
#define FFT8192_FREQ_GRID   500.0f      //频率栅格(500Hz整数倍)

// 外部全局变量声明
extern fft32_t fft8192_real[];      //实部缓冲
extern fft32_t fft8192_imag[];      //虚部缓冲
extern fft32_t fft8192_mag[];       //幅度谱
extern fft32_t fft8192_window[];    //汉宁窗系数
extern fft32_t fft8192_twiddle_cos[]; //预计算cos旋转因子
extern fft32_t fft8192_twiddle_sin[]; //预计算sin旋转因子
extern fft32_t fft8192_phase[];      //相位谱
extern fft32_t fft8192_magcorr[];     //补偿后幅度谱(每个bin)
extern FFT8192_Peak_t fft8192_peaks[]; //峰值数组
extern uint16_t fft8192_peakCnt;        //峰值数量
extern fft32_t fft8192_vrms;            //真有效值(mV)
extern fft32_t fft8192_vpp;             //峰峰值(mV)
extern fft32_t fft8192_out_freq[3];     //输出频率数组(Hz)
extern fft32_t fft8192_out_mag[3];      //输出幅度数组(mV)
extern fft32_t fft8192_out_phase[3];    //输出相位数组(rad)

// 函数声明
void FFT8192_Init(void);
void FFT8192_GenerateWindow(void);
void FFT8192_Process(const uint16_t* adc_in);
void FFT8192_FindPeaks(void);
void FFT8192_HarmonicMatch(void);
fft32_t FFT8192_PhaseDelay(fft32_t freq);   //滤波器相位延迟(弧度)
fft32_t FFT8192_FilterGain(fft32_t freq);   //滤波器增益补偿(1/gain)

#endif
