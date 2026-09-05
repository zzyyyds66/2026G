
#include "arm_math.h"
#include "math.h"
#define MAX_GROUP_COMPONENT 3U  //最多基波+2谐波

// 峰值结构体，与fft.c保持一致
typedef struct
{
    float32_t freq;
    float32_t mag;
		float32_t phase;
    uint8_t used;        //新增：是否被谐波组占用
} PeakInfo_t;
typedef struct
{
    uint8_t cnt;                                   //组内谱线数量 1~3
    PeakInfo_t peaks[MAX_GROUP_COMPONENT];
    float32_t sumMag;                              //幅值总和，择优使用
} HarmonicGroup_t;

typedef struct
{
    float32_t rms_normalized;       /* 归一化有效值（交流） */
    float32_t maximum_normalized;   /* 归一化峰值 */
    float32_t minimum_normalized;   /* 归一化谷值 */
    float32_t vpp_normalized;       /* 归一化峰峰值 */
    uint16_t  valid_input_samples;  /* 有效输入点数（不含边界丢弃） */
    uint16_t  interpolated_samples; /* 插值后总点数 */
} SignalMeasurement_t;

//==================== 外部全局变量声明 ====================
extern arm_cfft_instance_f32 scfft;
extern float32_t INPUT[];
extern float32_t FFT_INPUT[];
extern float32_t FFT_OUTPUT[];
extern float32_t magSpectrum[];
extern float32_t window_buffer[];
extern uint16_t adc_buf[];
extern float32_t adc_output[];

extern PeakInfo_t peakBuf[];
extern uint16_t peakCnt;
extern uint8_t peakOverflowFlag;

//==================== 函数声明 ====================
void RemoveDC(float32_t *INPUT);
void ADC_ArrayToVoltage(uint16_t *adc_in, float32_t *adc_out, uint16_t len);
void Get_Data(float32_t *wave, uint16_t len, float32_t *Vpp, float32_t *Vrms);
void FFT_Init();
void Generate_hanning_window();
void FFT_Process(float32_t *fft_in);
void SortPeakByFreq(PeakInfo_t *p, uint16_t num);
float32_t ParabolicInterp(float32_t y0, float32_t y1, float32_t y2, float32_t *delta_out);
uint32_t ClusterAndRoundFrequency(PeakInfo_t *inPeak, uint16_t inCnt, PeakInfo_t *outPeak);
uint8_t Signal_Process(float32_t *Vpp, float32_t *Vrms, float32_t *fA, float32_t *FREQ, float32_t *Mag);
uint16_t Get_Data_Average(float32_t *wave, uint16_t len, float32_t *Vpp_out, float32_t *Vrms_out);
void HarmonicMatch(PeakInfo_t *peakBuf, uint16_t peakCnt, HarmonicGroup_t *bestGroup);
uint8_t SignalMeasure_TimeDomain(const float32_t *input, uint16_t length, SignalMeasurement_t *result);