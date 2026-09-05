#ifndef ADC_DAC_H
#define ADC_DAC_H

#include "main.h"

//采样点数
#define SAMPLE_POINTS 8192 
//滤除点数
#define ADC_DISCARD_POINTS 16

// 全局变量
extern volatile uint8_t adc_done_flag;
extern uint16_t adc_buffer[SAMPLE_POINTS + ADC_DISCARD_POINTS];
extern uint16_t dac_buffer[SAMPLE_POINTS];
extern uint16_t adc_data[SAMPLE_POINTS];

// 接口函数
void ADCDAC_Init(void);
void ADCDAC_StartADC(void);
void ADCDAC_StopADC(void);
void ADCDAC_UpdateDAC(void);
void ADCDAC_ADCData(void);
void ADCDAC_ADCCopyDAC(void);


#endif /* ADC_DAC_H */
