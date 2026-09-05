#include "main.h"
#include "adc.h"
#include "tim.h"
#include "adcdac.h"
#include <string.h>

// 全局变量
volatile uint8_t adc_done_flag = 0;
uint16_t adc_buffer[SAMPLE_POINTS + ADC_DISCARD_POINTS] __attribute__((aligned(4))); //原始ADC采样数据
uint16_t adc_data[SAMPLE_POINTS] __attribute__((aligned(4)));    //处理后的ADC数据
uint16_t dac_buffer[SAMPLE_POINTS] __attribute__((aligned(4)));  //DAC输出缓冲区

// 初始化（只校准一次）
void ADCDAC_Init(void) {
    // ADC校准（只做一次）
    HAL_ADCEx_Calibration_Start(&hadc1, ADC_CALIB_OFFSET, ADC_SINGLE_ENDED);
    HAL_ADCEx_Calibration_Start(&hadc1, ADC_CALIB_OFFSET_LINEARITY, ADC_SINGLE_ENDED);
    // 清空缓冲
    memset(adc_buffer, 0, sizeof(adc_buffer));
    memset(dac_buffer, 0, sizeof(dac_buffer));
    memset(adc_data, 0, sizeof(adc_data));
    // 启动TIM6
    HAL_TIM_Base_Start(&htim6);
}

// 启动ADC阻塞采样采满缓冲区
void ADCDAC_StartADC(void) {
    adc_done_flag = 0;
    HAL_ADC_Start_DMA(&hadc1, (uint32_t*)adc_buffer, SAMPLE_POINTS + ADC_DISCARD_POINTS);
    while (!adc_done_flag) HAL_Delay(1);
}

// 停止ADC采样
void ADCDAC_StopADC(void) {
    HAL_ADC_Stop_DMA(&hadc1);
//    ADC可关可不关吧可能？
//    HAL_ADC_Stop(&hadc1);
}

// 丢弃不稳定点，得到adcdata数据
void ADCDAC_ADCData(void) {
    //截取有效数据存入处理缓冲区
    memcpy(adc_data, adc_buffer + ADC_DISCARD_POINTS, SAMPLE_POINTS * sizeof(uint16_t));
}

void ADCDAC_ADCCopyDAC(void) {
    //直通下同步数据到DAC输出缓存
    memcpy(dac_buffer, adc_data, SAMPLE_POINTS * sizeof(uint16_t));
}

/*
// 更新DAC输出
void ADCDAC_UpdateDAC(void) {
    HAL_DAC_Stop_DMA(&hdac1, DAC_CHANNEL_1);
    HAL_DAC_Start_DMA(&hdac1, DAC_CHANNEL_1, (uint32_t*)dac_buffer, SAMPLE_POINTS, DAC_ALIGN_12B_R);
}
*/
// ADC DMA完成回调
void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef* hadc) {
    if (hadc->Instance == ADC1) {
        adc_done_flag = 1;
    }
}

// DAC DMA完成回调（循环模式下不会触发）
//void HAL_DAC_ConvCpltCallbackCh1(DAC_HandleTypeDef* hdac) {}
