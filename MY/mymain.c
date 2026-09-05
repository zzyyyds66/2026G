#include "main.h"
#include "mymain.h"
#include "adcdac.h"
#include "arm_math.h"
#include "fft.h"

// 工作模式
WorkMode mode = MODE_IDLE;
float Vpp = 0;
float Vrms = 0;
float fA = 0;

// 主循环处理
void mymain_process(void) {

    switch (mode) {
    case MODE_IDLE:
        break;

    case MODE_PASS://直通模式，但是为了简单直接使用了单次采集，符合后续示波器的完成
        ADCDAC_StartADC();
        ADCDAC_StopADC();
        ADCDAC_ADCData();
//        ADCDAC_ADCCopyDAC();
//        ADCDAC_UpdateDAC();
        break;

    case MODE_SCOP:
        // TODO: 加窗 + FFT + 参数计算 + 串口发送
        ADCDAC_StartADC();
        ADCDAC_StopADC();
        ADCDAC_ADCData();
//		    ADCDAC_ADCCopyDAC();
//		    ADCDAC_UpdateDAC();
		
		    memcpy(adc_buf, adc_data, SAMPLE_POINTS * sizeof(uint16_t)); 	
		
        break;
    
    default:
        break;
    }
}