#ifndef DWIN_H
#define DWIN_H

#include "stm32h7xx_hal.h"

#ifdef __cplusplus
extern "C" {
#endif

    // 取最大值
#define MAX(a, b)      ((a) > (b) ? (a) : (b))
// 取最小值
#define MIN(a, b)      ((a) < (b) ? (a) : (b))
// 最小最大限幅
#define MINMAX(x, low, high)   MAX((low), MIN((x), (high)))

/* 按键类型 */
    typedef enum {
        DWIN_KEY_NONE = 0,
        DWIN_KEY_1 = 1,
        DWIN_KEY_2 = 2
    } DWIN_KeyTypeDef;

    /* 全局变量 */
    extern volatile DWIN_KeyTypeDef dwin_key;
    extern uint8_t T_mode;      // 周期模式(0:单周期, 1:三周期)
    extern uint8_t F_mode;      // 时频域模式(0:时域, 1:频域)
    extern uint8_t wave_src;    // 波形源(0:基础合成, 1:直通ADC, 2:带相位合成, 3:LS拟合)
    extern float dw_phase[3];   // 各次谐波相位(rad)
    //extern volatile uint8_t dwin_reply_timeout;

#define TEXT1_ADDR          0x5300      // 文本框VP地址
#define TEXT2_ADDR          0x5400      // 文本框VP地址
#define TEXT3_ADDR          0x5500      // 文本框VP地址
#define TEXT_MAX_LEN        100         // 文本框最大长度(字节)

/* 公共接口 */
    HAL_StatusTypeDef DWIN_Init(UART_HandleTypeDef* huart);
    void DWIN_ParseFrame(uint16_t size);
    void DWIN_ClearText(uint16_t ch);
    void DWIN_SendText(uint16_t ch, const char* fmt, ...);
    void DWIN_ClearWave(void);
    void DWIN_SendWave(uint8_t channel, const uint16_t* data, uint16_t num);
    void DWIN_Process(void);
    void DWIN_Main(void);

#ifdef __cplusplus
}
#endif

#endif