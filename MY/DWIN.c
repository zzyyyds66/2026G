#include "dwin.h"
#include "fft.h"
#include "fft8192.h"
#include "adcdac.h"
#include "arm_math.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>

extern UART_HandleTypeDef huart1;

#define M_PI 3.14159265359f

/* 协议常量 */
#define CMD_WRITE           0x82        // 写命令
#define CMD_READ            0x83        // 读应答(按键)

/* 地址配置(来自地址及键值.txt) */
#define BTN_ADDR_H          0x60        // 按键VP地址高字节
#define BTN_ADDR_L          0x00        // 按键VP地址低字节
#define BTN_WORD_CNT        0x01        // 按键字计数
#define BTN_DATA_H          0x00        // 按键数据高字节


/* 缓冲区配置 */
#define RX_BUF_SIZE         9           /* 串口接收缓冲区长度，同时刚好是按键返回帧的长度 */
#define TX_BUF_SIZE         256         /* 串口发送缓冲区长度 */
#define TX_WAVE_MAX_POINTS  64          /* 单次串口发送最大点数 */
#define WAVE_POINTS         768         /* 显示波形数据横轴点数 */
#define VPP_POINTS          4096        /* Vpp高分辨率采样点数 */
#define WAVE_AMPL           2000        /* 显示波形数据纵轴幅度 */
#define ADC_SAMPLE_RATE     3000000.0f  /* ADC采样频率: 3MHz */


/* 全局变量 */
volatile DWIN_KeyTypeDef dwin_key = DWIN_KEY_1;
//volatile uint8_t dwin_reply_timeout = 0;
uint16_t wave_data[WAVE_POINTS];  // 波形数据缓冲区(任意长度, 自动分块发送)
uint8_t T_mode = 0; // 显示时域周期数(0:单周期, 1:三周期)
uint8_t F_mode = 0; // 显示时频域模式(0:时域, 1:频域)
uint8_t wave_src = 0; // 波形源(0:直通ADC, 1:数字合成)
uint8_t start_run = 0; // 开始更新波形

float dw_umax = -1.0f; // 最大值
float dw_umin = 1.0f; // 最小值
float dw_upp = 0.0f;  // 峰峰值
float dw_urms = 0.0f; // 有效值
float dw_f1 = 1000.0f;   // 基频率
float dw_f2 = 3000.0f;   // 谐波1频率
float dw_f3 = 5000.0f;   // 谐波2频率
float dw_u1 = 100.0f;   // 基幅度
float dw_u2 = 33.0f;   // 谐波1幅度
float dw_u3 = 20.0f;   // 谐波2幅度
float dw_p1 = 0.0f;   // 基相位
float dw_p2 = 0.0f;   // 谐波1相位
float dw_p3 = 0.0f;   // 谐波2相位

float dw_sin[WAVE_POINTS] = { 0 };   // 标准正弦波数据[-1,1]
uint16_t k2 = 0.0f; // 谐波1相对基波的频率比
uint16_t k3 = 0.0f; // 谐波2相对基波的频率比
float dw_phase[3] = { 0 }; // 各次谐波相位(rad)
float sum_amp = 0.0f; // 叠加总峰值(mV)
float wave1[WAVE_POINTS] = { 0 };   // 基波数据[-1,1]
float wave2[WAVE_POINTS] = { 0 };   // 谐波1数据[-1,1]
float wave3[WAVE_POINTS] = { 0 };   // 谐波2数据[-1,1]
float waveall[WAVE_POINTS] = { 0 };   // 总波形数据[-1,1]
uint8_t numx = 15;// 频谱谱线数
uint16_t deltax = 15;// 频谱谱线间隙

/* 调试用全局变量 */
float freq_buf[10][3] = { 0 };   // 10次采样频点记录
float mag_buf[10][3] = { 0 };    // 10次采样幅值记录
float vpp_buf[10] = { 0 };       // 10次时域Vpp记录(mV)
float vrms_buf[10] = { 0 };      // 10次时域Vrms记录(mV)
float x_buf[10][6] = { 0 };      // 10次LS系数记录
float x_mean[6] = { 0 };         // LS系数trimmed mean
float dw_upp_td = 0.0f;          // 时域Vpp trimmed mean(mV)
float dw_urms_td = 0.0f;         // 时域Vrms trimmed mean(mV)
uint8_t valid_cnt = 0;            // 有效检测次数
uint8_t vote_cnt[10] = { 0 };     // 基频投票计数
float mf1 = 0.0f;                 // 主频(投票胜出)
float wave_fit[VPP_POINTS] = { 0 };   // 高分辨率单周期波形(归一化[-1,1])
float fit_upp = 0.0f;              // LS拟合幅度Vpp(mV), 供对比
float fit_amp[3] = { 0.0f };
float fit_phase[3] = { 0.0f };

float abc = 0.98f;

/* 内部变量 */
static UART_HandleTypeDef* g_huart;
static uint8_t rx_buffer[RX_BUF_SIZE];
static uint8_t tx_buffer[TX_BUF_SIZE];
static char fmt_buffer[TEXT_MAX_LEN + 1];  // 格式化临时缓冲区
//static volatile uint8_t waiting_reply = 0;
//static volatile uint32_t send_tick = 0;


//初始化迪文屏驱动
HAL_StatusTypeDef DWIN_Init(UART_HandleTypeDef* huart) {
    if (huart == NULL) return HAL_ERROR;
    g_huart = huart;
    HAL_UARTEx_ReceiveToIdle_IT(huart, rx_buffer, RX_BUF_SIZE);

    // 发送复位帧: 5a a5 07 82 00 04 55 aa 5a a5
    uint8_t reset_frame[] = { 0x5A, 0xA5, 0x07, 0x82, 0x00, 0x04, 0x55, 0xAA, 0x5A, 0xA5 };
    HAL_UART_Transmit(huart, reset_frame, sizeof(reset_frame), 1000);
    HAL_Delay(10);

    // 初始化正弦波数据
    for (int i = 0; i < WAVE_POINTS; i++) {
        dw_sin[i] = sinf(i * 2 * M_PI / WAVE_POINTS);
    }

    // 清除初始化期间可能产生的错误
    //__HAL_UART_CLEAR_OREFLAG(huart);
    //__HAL_UART_CLEAR_FEFLAG(huart);

        //开启空闲串口中断
        //return HAL_UARTEx_ReceiveToIdle_IT(huart, rx_buffer, RX_BUF_SIZE);
    return 1;
}

/**
 * @brief  解析接收帧(在串口空闲中断回调中调用)
 * @param  size: 实际接收字节数
 * @note   按键帧固定9字节: 5A A5 06 83 60 00 01 00 XX
 */
void DWIN_ParseFrame(uint16_t size) {
    // 严格检查: 只处理9字节的按键帧
    if (size != RX_BUF_SIZE) goto restart;

    // 校验帧头、长度、命令字、地址、字计数、数据高字节
    if (rx_buffer[0] != 0x5A || rx_buffer[1] != 0xA5 ||
        rx_buffer[2] != 0x06 || rx_buffer[3] != CMD_READ ||
        rx_buffer[4] != BTN_ADDR_H || rx_buffer[5] != BTN_ADDR_L ||
        rx_buffer[6] != BTN_WORD_CNT || rx_buffer[7] != BTN_DATA_H) {
        goto restart;
    }

    // 提取键值
    uint8_t key_val = rx_buffer[8];
    if (key_val >= 1 && key_val <= 8) {
        dwin_key = (DWIN_KeyTypeDef)key_val;
    }

restart:
    HAL_UARTEx_ReceiveToIdle_IT(g_huart, rx_buffer, RX_BUF_SIZE);
}

//清空文本框
void DWIN_ClearText(uint16_t ch) {
    // 5A A5 05 82 53/54/55 00 00 00
    uint8_t chhigh = ch >> 8;
    uint8_t frame[] = { 0x5A, 0xA5, 0x05, 0x82, chhigh, 0x00, 0x00, 0x00 };
    HAL_UART_Transmit(g_huart, frame, sizeof(frame), 100);
    HAL_Delay(10);
}


//发送文本(格式化输入: 支持%d, %f等占位符)
void DWIN_SendText(uint16_t ch, const char* fmt, ...) {
    if (fmt == NULL) return;

    // 格式化字符串
    va_list args;
    va_start(args, fmt);
    vsnprintf(fmt_buffer, sizeof(fmt_buffer), fmt, args);
    va_end(args);

    // 清空文本框
    DWIN_ClearText(ch);
    HAL_Delay(5);

    // 计算文本长度
    uint16_t len = 0;
    while (fmt_buffer[len] && len < TEXT_MAX_LEN) len++;

    // 构建写入帧: 5A A5 [len] 82 [addr] [word_cnt] [data...]
    uint16_t data_len = 3 + TEXT_MAX_LEN;  // addr(2) + word_cnt(1) + data(100)
    uint16_t total_len = 4 + data_len;
    if (total_len > TX_BUF_SIZE) return;

    tx_buffer[0] = 0x5A;
    tx_buffer[1] = 0xA5;
    tx_buffer[2] = data_len;
    tx_buffer[3] = 0x82;
    tx_buffer[4] = (uint8_t)(ch >> 8);
    tx_buffer[5] = (uint8_t)(ch & 0xFF);

    for (uint16_t i = 0; i < TEXT_MAX_LEN; i++) {
        tx_buffer[6 + i] = (i < len) ? (uint8_t)fmt_buffer[i] : 0x00;
    }

    HAL_UART_Transmit(g_huart, tx_buffer, total_len, 1000);
}

//清空波形通道
void DWIN_ClearWave(void) {
    // 5A A5 05 82 03 08 00 00
//    uint8_t frame[] = { 0x5A, 0xA5, 0x05, 0x82, 0x03, 0x08, 0x00, 0x00 };
//		uint8_t frame[] = { 0x5A, 0xA5, 0x05, 0x82, 0x04, 0x00, 0x00, 0x00 };
    uint8_t frame[] = { 0x5A, 0xA5, 0x05, 0x82, 0x03, 0x09, 0x00, 0x00 };
    HAL_UART_Transmit(g_huart, frame, sizeof(frame), 100);
    HAL_Delay(10);
}

/**
 * @brief  发送波形数据(清空+分块填充)
 * @param  channel: 波形通道号(0x00-0x07)
 * @param  data: 波形数据数组
 * @param  num: 数据点数量(任意长度, 自动分块)
 */
void DWIN_SendWave(uint8_t channel, const uint16_t* data, uint16_t num) {
    if (data == NULL || num == 0) return;

    // 清空通道
    DWIN_ClearWave();

    // 分块发送, 每块最多TX_WAVE_MAX_POINTS个点
    uint16_t offset = 0;
    while (offset < num) {
        //本轮发送16位数据的点数
        uint16_t chunk = (num - offset > TX_WAVE_MAX_POINTS) ? TX_WAVE_MAX_POINTS : (num - offset);

        // 帧1: 启动曲线缓冲区写 5A A5 03 82 00 10
        //uint8_t cmd_start[] = { 0x5A, 0xA5, 0x03, 0x82, 0x00, 0x10 };
        //HAL_UART_Transmit(g_huart, cmd_start, sizeof(cmd_start), 100);

        // 帧: 5A A5 [len] 82 03 10 5A A5 01 00 [channel] [points] [data...]
        uint16_t data_len = 9 + chunk * 2;
        uint16_t total_len = 3 + data_len;
        if (total_len > TX_BUF_SIZE) break;

        tx_buffer[0] = 0x5A;
        tx_buffer[1] = 0xA5;
        tx_buffer[2] = data_len;
        tx_buffer[3] = 0x82;
        tx_buffer[4] = 0x03;
        tx_buffer[5] = 0x10;
        tx_buffer[6] = 0x5A;
        tx_buffer[7] = 0xA5;
        tx_buffer[8] = 0x01;
        tx_buffer[9] = 0x00;
        tx_buffer[10] = channel;
        tx_buffer[11] = (uint8_t)chunk;


        for (uint16_t i = 0, buf_pos = 12; i < chunk; i++) {
            uint16_t val = data[offset + i];
            tx_buffer[buf_pos++] = val >> 8;
            tx_buffer[buf_pos++] = val & 0xFF;
        }

        HAL_UART_Transmit(g_huart, tx_buffer, total_len, 1000);
        offset += chunk;
        HAL_Delay(20);
    }

    //    waiting_reply = 1;
    //    send_tick = HAL_GetTick();
    //    dwin_reply_timeout = 0;
}

/**
 * @brief  检查回复超时(在主循环中调用)

void DWIN_Process(void) {
    if (waiting_reply) {
        if ((HAL_GetTick() - send_tick) > REPLY_TIMEOUT_MS) {
            dwin_reply_timeout = 1;
            waiting_reply = 0;
        }
    }
}
 */


 //串口空闲中断回调
void HAL_UARTEx_RxEventCallback(UART_HandleTypeDef* huart, uint16_t Size) {
    if (huart->Instance == USART1) {
        DWIN_ParseFrame(Size);
    }
}

//串口错误中断回调
void HAL_UART_ErrorCallback(UART_HandleTypeDef* huart) {
    if (huart->Instance == USART1) {
        // 清空溢出、帧错误标记
        __HAL_UART_CLEAR_OREFLAG(huart);
        __HAL_UART_CLEAR_FEFLAG(huart);
        // 重新开启空闲接收
        HAL_UARTEx_ReceiveToIdle_IT(huart, rx_buffer, RX_BUF_SIZE);
    }
}

/**
 * @brief  直通模式: ADC数据直接映射到显示波形, 自动缩放上/过零点对齐
 * @note   使用FFT检测的周期计算显示窗口, 寻找第一个上升沿过零点作为起点
 */
 /*
 static void DWIN_UpdateWaveDirect(void) {
     // 1. 由FFT频率计算周期(采样点数)
     uint16_t period_samples;
     if (dw_f1 > 0.01f && dw_f1 < ADC_SAMPLE_RATE / 2.0f) {
         period_samples = (uint16_t)(ADC_SAMPLE_RATE / dw_f1);
         period_samples = MAX(period_samples, 8);
     }
     else {
         period_samples = 1000;
     }

     // 2. 根据T_mode确定显示周期数: 0→单周期, 1→三周期
     uint8_t cycles = T_mode ? 3 : 1;
     uint32_t display_samples = (uint32_t)period_samples * cycles;
     display_samples = MINMAX(display_samples, WAVE_POINTS, SAMPLE_POINTS);

     // 3. 在有效区间内搜索信号范围, 用于自动缩放
     uint16_t sig_min = 0xFFFF, sig_max = 0;
     uint16_t search_limit = MIN(display_samples, SAMPLE_POINTS);

     for (uint16_t i = 0; i < search_limit; i++) {
         uint16_t val = adc_data[i];
         sig_min = MIN(sig_min, val);
         sig_max = MAX(sig_max, val);
     }

     // 4. 计算信号中心和半幅值范围
     int32_t sig_center = ((int32_t)sig_min + (int32_t)sig_max) / 2;
     int32_t sig_half = ((int32_t)sig_max - (int32_t)sig_min + 1) / 2;
     if (sig_half == 0) sig_half = 1;

     // 5. 搜索第一个上升沿过零点 (adc_data跨过sig_center的位置)
     uint16_t zero_cross = 0;
     for (uint16_t i = 1; i < search_limit; i++) {
         int32_t prev = (int32_t)adc_data[i - 1];
         int32_t curr = (int32_t)adc_data[i];
         if (prev <= sig_center && curr > sig_center) {
             int32_t delta = curr - prev;
             if (delta > 0) {
                 int32_t frac = (sig_center - prev) * 1000 / delta;
                 zero_cross = i - 1 + (uint16_t)(frac / 1000);
             }
             break;
         }
     }
     if (zero_cross >= display_samples) zero_cross = 0;

     // 6. 从过零点开始降采样, 映射到显示范围[0, 2*WAVE_AMPL]
     for (uint16_t i = 0; i < WAVE_POINTS; i++) {
         uint32_t src_idx = zero_cross + (uint32_t)i * display_samples / WAVE_POINTS;
         if (src_idx >= SAMPLE_POINTS) src_idx = SAMPLE_POINTS - 1;
         int32_t centered = (int32_t)adc_data[src_idx] - sig_center;
         int32_t display = centered * (int32_t)WAVE_AMPL / sig_half + (int32_t)WAVE_AMPL;
         wave_data[i] = (uint16_t)MINMAX(display, 0, 2 * WAVE_AMPL);
     }

     // 7. 更新测量值(转mV, 复用FFT8192已算的Vrms/Vpp)
     dw_umin = (float)sig_min * 3300.0f / 65535.0f;
     dw_umax = (float)sig_max * 3300.0f / 65535.0f;
     dw_upp = fft8192_vpp;
     dw_urms = fft8192_vrms;
 }
 */

 /**
  * @brief  带相位的完善数字合成: 使用FFT检测到的频率、幅度、相位进行波形重构
  * @note   利用精确的频率比(浮点数)和FFT相位, 合成更接近原始信号的波形
  */
  // 三周期降采样: 将单周期波形压缩为三周期显示
static void ApplyTMode(void) {
    for (uint16_t i = 0; i < WAVE_POINTS / 3; i++)
        wave_data[i] = wave_data[3 * i];
    for (uint16_t i = WAVE_POINTS / 3; i < WAVE_POINTS * 2 / 3; i++)
        wave_data[i] = wave_data[i - WAVE_POINTS / 3];
    for (uint16_t i = WAVE_POINTS * 2 / 3; i < WAVE_POINTS; i++)
        wave_data[i] = wave_data[i - WAVE_POINTS * 2 / 3];
}

// 基础数字合成: 基波+谐波叠加(无相位), 使用dw_sin查表
/*
static void DWIN_UpdateWaveBasic(void) {
    dw_umin = 999.0f;
    dw_umax = -999.0f;
    dw_urms = 0.0f;

    for (uint16_t i = 0; i < WAVE_POINTS; i++) {
        wave1[i] = dw_u1 * dw_sin[i];
        wave2[i] = dw_u2 * dw_sin[(i * k2) % WAVE_POINTS];
        wave3[i] = dw_u3 * dw_sin[(i * k3) % WAVE_POINTS];
        waveall[i] = (wave1[i] + wave2[i] + wave3[i]) / sum_amp;
        dw_umin = MIN(dw_umin, waveall[i]);
        dw_umax = MAX(dw_umax, waveall[i]);
        dw_urms += waveall[i] * waveall[i];
        wave_data[i] = WAVE_AMPL * (waveall[i] + 1);
    }
    dw_upp = (dw_umax - dw_umin) * sum_amp;
    dw_urms = sqrtf(dw_urms / WAVE_POINTS) * sum_amp;

    if (T_mode) ApplyTMode();
}
*/

// 带相位数字合成: 基波+谐波叠加(含相位), 使用dw_sin查表+相位偏移
/*
static void DWIN_UpdateWaveSynth(void) {
    // 相位转索引偏移: shift = phase * N / (2π)
    int16_t shift1 = (int16_t)((dw_p2 - dw_p1) * WAVE_POINTS / (2.0f * PI));
    int16_t shift2 = (int16_t)((dw_p3 - dw_p1) * WAVE_POINTS / (2.0f * PI));

    dw_umin = 999.0f;
    dw_umax = -999.0f;
    dw_urms = 0.0f;

    for (uint16_t i = 0; i < WAVE_POINTS; i++) {
        int32_t idx1 = (int32_t)(i * k2) + shift1;
        int32_t idx2 = (int32_t)(i * k3) + shift2;
        wave1[i] = dw_u1 * dw_sin[i];
        wave2[i] = dw_u2 * dw_sin[((idx1 % WAVE_POINTS) + WAVE_POINTS) % WAVE_POINTS];
        wave3[i] = dw_u3 * dw_sin[((idx2 % WAVE_POINTS) + WAVE_POINTS) % WAVE_POINTS];
        waveall[i] = (wave1[i] + wave2[i] + wave3[i]) / sum_amp;
        dw_umin = MIN(dw_umin, waveall[i]);
        dw_umax = MAX(dw_umax, waveall[i]);
        dw_urms += waveall[i] * waveall[i];
        wave_data[i] = WAVE_AMPL * (waveall[i] + 1);
    }
    dw_upp = (dw_umax - dw_umin) * sum_amp;
    dw_urms = sqrtf(dw_urms / WAVE_POINTS) * sum_amp;

    if (T_mode) ApplyTMode();
}
*/

// 牛顿迭代精估极值: 在t0附近求s'(t)=0, 返回精估后的s(t)
static float NewtonPeak(float t0, const float* amps, const float* freqs, const float* phases, uint8_t n) {
    float t = t0;
    for (uint8_t iter = 0; iter < 3; iter++) {
        float sp = 0, spp = 0;
        for (uint8_t j = 0; j < n; j++) {
            float w = 2.0f * PI * freqs[j];
            float wt = w * t + phases[j];
            sp += amps[j] * w * cosf(wt);
            spp -= amps[j] * w * w * sinf(wt);
        }
        if (fabsf(spp) < 1e-12f) break;   // 防除零
        t -= sp / spp;
    }
    float s = 0;
    for (uint8_t j = 0; j < n; j++)
        s += amps[j] * sinf(2.0f * PI * freqs[j] * t + phases[j]);
    return s;
}

// 最小二乘拟合波形: 使用10次采样trimmed mean后的LS系数合成
static void DWIN_UpdateWaveFit(void) {
    if (dw_f1 < 0.01f) return;   // 无有效信号, 跳过拟合
    uint8_t n_freq = (dw_f2 > 0.01f) ? ((dw_f3 > 0.01f) ? 3 : 2) : 1;
    float freqs[3] = { dw_f1, dw_f2, dw_f3 };

    // 从x_mean提取幅度和相位供调试(已在采样循环内补偿滤波器相位)
    memset(fit_amp, 0, sizeof(fit_amp));
    memset(fit_phase, 0, sizeof(fit_phase));
    for (uint8_t j = 0; j < n_freq; j++) {
        fit_amp[j] = sqrtf(x_mean[2 * j] * x_mean[2 * j] + x_mean[2 * j + 1] * x_mean[2 * j + 1]) * 3300.0f;
        fit_phase[j] = atan2f(x_mean[2 * j + 1], x_mean[2 * j]);
    }

    // 高分辨率波形+Vpp: 直接用x_mean合成(已补偿+对齐)
    {
        float s_max = -1e9f, s_min = 1e9f;
        for (uint16_t i = 0; i < VPP_POINTS; i++) {
            float t = (float)i / VPP_POINTS / dw_f1;
            float val = 0.0f;
            for (uint8_t j = 0; j < n_freq; j++) {
                float wt = 2.0f * PI * freqs[j] * t;
                val += x_mean[2 * j] * sinf(wt) + x_mean[2 * j + 1] * cosf(wt);
            }
            val *= 3300.0f;
            if (val > s_max) s_max = val;
            if (val < s_min) s_min = val;
            wave_fit[i] = val / sum_amp;
        }
        dw_upp = s_max - s_min;
        fit_upp = dw_upp;
    }

    // 7. 抽768点显示
    for (uint16_t i = 0; i < WAVE_POINTS; i++)
        wave_data[i] = WAVE_AMPL * (wave_fit[(uint32_t)i * VPP_POINTS / WAVE_POINTS] + 1.0f);

    // 8. 三周期模式: 复用ApplyTMode
    if (T_mode) ApplyTMode();

    // 9. 有效值: FFT幅度解析式 Σ(U²/2) 开方
    dw_urms = sqrtf((dw_u1 * dw_u1 + dw_u2 * dw_u2 + dw_u3 * dw_u3) * 0.5f);

}

// 截尾平均: 去最大最小后取平均, 不足3个直接取平均, 0个返回0
static float TrimmedMean(const float* vals, uint8_t cnt) {
    if (cnt == 0) return 0.0f;
    if (cnt < 3) {
        float sum = 0;
        for (uint8_t i = 0; i < cnt; i++) sum += vals[i];
        return sum / cnt;
    }
    uint8_t imax = 0, imin = 0;
    for (uint8_t i = 1; i < cnt; i++) {
        if (vals[i] > vals[imax]) imax = i;
        if (vals[i] < vals[imin]) imin = i;
    }
    float sum = 0; uint8_t rem = 0;
    for (uint8_t i = 0; i < cnt; i++) {
        if (i == imax || i == imin) continue;
        sum += vals[i]; rem++;
    }
    return sum / rem;
}

void DWIN_Main(void) {
    // 按键处理
    if (dwin_key != DWIN_KEY_NONE) {
        start_run = 1;
        switch (dwin_key) {
        case DWIN_KEY_1:
            if (F_mode == 0) {
                // 时域模式: 切换周期数
                T_mode = 1 - T_mode;
            }
            //            else {
                            // 频域模式: 循环切换波形源 (基础→直通→带相位→LS拟合)
            //                wave_src = (wave_src + 1) % 4;
            //            }
            break;
        case DWIN_KEY_2:
            //切换时域频域
            F_mode = 1 - F_mode;
        default:
            HAL_Delay(10);
            break;
        }
        dwin_key = DWIN_KEY_NONE;
    }
    else {
        HAL_Delay(100);
    }

    if (start_run == 1) {
        start_run = 0;

        // 10次采样+FFT, 记录频点幅值组合, 投票去偶发误差
        valid_cnt = 0;
        for (uint8_t n = 0; n < 10; n++) {
            ADCDAC_StartADC(); //阻塞式采样
            ADCDAC_StopADC();  //停止采样
            ADCDAC_ADCData();
            FFT8192_Process(adc_data);
            if (fft8192_out_freq[0] > 0.01f) {   // 过滤无效检测
                for (uint8_t j = 0; j < 3; j++) {
                    freq_buf[valid_cnt][j] = fft8192_out_freq[j];
                    mag_buf[valid_cnt][j] = fft8192_out_mag[j];
                }
                vpp_buf[valid_cnt] = fft8192_vpp;
                vrms_buf[valid_cnt] = fft8192_vrms;
                // LS拟合系数记录: 用当前FFT频率对当前ADC数据做LS
                {
                    float lf[3] = { fft8192_out_freq[0], fft8192_out_freq[1], fft8192_out_freq[2] };
                    uint8_t l_nfreq = (lf[0] > 0.01f) ? ((lf[1] > 0.01f) ? ((lf[2] > 0.01f) ? 3 : 2) : 1) : 0;
                    if (l_nfreq > 0) {
                        float l_mean = 0.0f;
                        for (uint16_t i = 0; i < SAMPLE_POINTS; i++)
                            l_mean += ((float)adc_data[i] / 32767.5f - 1.0f - l_mean) / (float)(i + 1);
                        float l_A[6][6] = { 0 }, l_b[6] = { 0 };
                        uint8_t l_nparam = 2 * l_nfreq;
                        for (uint16_t i = 0; i < SAMPLE_POINTS; i++) {
                            float sig = (float)adc_data[i] / 32767.5f - 1.0f - l_mean;
                            float t = (float)i / ADC_SAMPLE_RATE;
                            float basis[6];
                            for (uint8_t j = 0; j < l_nfreq; j++) {
                                float wt = 2.0f * PI * lf[j] * t;
                                basis[2 * j] = sinf(wt);
                                basis[2 * j + 1] = cosf(wt);
                            }
                            for (uint8_t r = 0; r < l_nparam; r++) {
                                l_b[r] += sig * basis[r];
                                for (uint8_t c = r; c < l_nparam; c++) l_A[r][c] += basis[r] * basis[c];
                            }
                        }
                        for (uint8_t r = 0; r < l_nparam; r++)
                            for (uint8_t c = 0; c < r; c++) l_A[r][c] = l_A[c][r];
                        float l_x[6] = { 0 };
                        for (uint8_t k = 0; k < l_nparam; k++) {
                            uint8_t piv = k;
                            for (uint8_t r = k + 1; r < l_nparam; r++)
                                if (fabsf(l_A[r][k]) > fabsf(l_A[piv][k])) piv = r;
                            if (piv != k) {
                                for (uint8_t c = 0; c < l_nparam; c++) { float tmp = l_A[k][c]; l_A[k][c] = l_A[piv][c]; l_A[piv][c] = tmp; }
                                float tmp = l_b[k]; l_b[k] = l_b[piv]; l_b[piv] = tmp;
                            }
                            for (uint8_t r = k + 1; r < l_nparam; r++) {
                                float f = l_A[r][k] / l_A[k][k];
                                for (uint8_t c = k; c < l_nparam; c++) l_A[r][c] -= f * l_A[k][c];
                                l_b[r] -= f * l_b[k];
                            }
                        }
                        for (int8_t r = l_nparam - 1; r >= 0; r--) {
                            l_x[r] = l_b[r];
                            for (uint8_t c = r + 1; c < l_nparam; c++) l_x[r] -= l_A[r][c] * l_x[c];
                            l_x[r] /= l_A[r][r];
                        }
                        // 1. 幅度补偿+相位补偿: 还原信号源幅度和相位
                        for (uint8_t j = 0; j < l_nfreq; j++) {
                            float gain_inv = FFT8192_FilterGain(lf[j]);   // 1/gain
                            l_x[2 * j] *= gain_inv;
                            l_x[2 * j + 1] *= gain_inv;
                            float delay = FFT8192_PhaseDelay(lf[j]);
                            float cs = cosf(delay), sn = sinf(delay);
                            float a = l_x[2 * j], b = l_x[2 * j + 1];
                            l_x[2 * j] = a * cs - b * sn;
                            l_x[2 * j + 1] = a * sn + b * cs;
                        }
                        // 2. 相位对齐: 以基波时间起点为参考, 各频率按次数旋转到φ0=0
                        float phi0 = atan2f(l_x[1], l_x[0]);   // 补偿后基波相位(纯信号源)
                        for (uint8_t j = 0; j < l_nfreq; j++) {
                            float angle = -phi0 * lf[j] / lf[0];
                            float cs = cosf(angle), sn = sinf(angle);
                            float a = l_x[2 * j], b = l_x[2 * j + 1];
                            x_buf[valid_cnt][2 * j] = a * cs - b * sn;
                            x_buf[valid_cnt][2 * j + 1] = a * sn + b * cs;
                        }
                    }
                }
                valid_cnt++;
            }
        }

        if (valid_cnt == 0) {
            dw_f1 = dw_f2 = dw_f3 = 0;
            dw_u1 = dw_u2 = dw_u3 = 0;
        }
        else {
            // 投票: 按完整频点组合(f1,f2,f3)分组, 滤除谐波匹配错误的采样
            uint8_t max_v = 0, major = 0;
            for (uint8_t n = 0; n < valid_cnt; n++) {
                vote_cnt[n] = 0;
                for (uint8_t m = 0; m < valid_cnt; m++)
                    if (fabsf(freq_buf[m][0] - freq_buf[n][0]) < 1.0f &&
                        fabsf(freq_buf[m][1] - freq_buf[n][1]) < 1.0f &&
                        fabsf(freq_buf[m][2] - freq_buf[n][2]) < 1.0f) vote_cnt[n]++;
                if (vote_cnt[n] > max_v) { max_v = vote_cnt[n]; major = n; }
            }
            mf1 = freq_buf[major][0];

            // 主组内 trimmed mean (主组: 频点组合一致, 1Hz容差)
            float vals[10];
            float* dw_u[3] = { &dw_u1, &dw_u2, &dw_u3 };
            for (uint8_t j = 0; j < 3; j++) {
                uint8_t cnt = 0;
                for (uint8_t n = 0; n < valid_cnt; n++)
                    if (fabsf(freq_buf[n][0] - freq_buf[major][0]) < 1.0f &&
                        fabsf(freq_buf[n][1] - freq_buf[major][1]) < 1.0f &&
                        fabsf(freq_buf[n][2] - freq_buf[major][2]) < 1.0f &&
                        mag_buf[n][j] > 0.01f) vals[cnt++] = mag_buf[n][j];
                *dw_u[j] = TrimmedMean(vals, cnt);
            }
            // 时域Vpp/Vrms同样做trimmed mean
            {
                uint8_t cnt = 0;
                for (uint8_t n = 0; n < valid_cnt; n++)
                    if (fabsf(freq_buf[n][0] - freq_buf[major][0]) < 1.0f &&
                        fabsf(freq_buf[n][1] - freq_buf[major][1]) < 1.0f &&
                        fabsf(freq_buf[n][2] - freq_buf[major][2]) < 1.0f &&
                        vpp_buf[n] > 0.01f) vals[cnt++] = vpp_buf[n];
                dw_upp_td = 0.5f * TrimmedMean(vals, cnt) * abc;
                cnt = 0;
                for (uint8_t n = 0; n < valid_cnt; n++)
                    if (fabsf(freq_buf[n][0] - freq_buf[major][0]) < 1.0f &&
                        fabsf(freq_buf[n][1] - freq_buf[major][1]) < 1.0f &&
                        fabsf(freq_buf[n][2] - freq_buf[major][2]) < 1.0f &&
                        vrms_buf[n] > 0.01f) vals[cnt++] = vrms_buf[n];
                dw_urms_td = TrimmedMean(vals, cnt);
            }
            // LS系数trimmed mean: 主组内对每个系数分别做trimmed mean
            for (uint8_t j = 0; j < 6; j++) {
                uint8_t cnt = 0;
                for (uint8_t n = 0; n < valid_cnt; n++)
                    if (fabsf(freq_buf[n][0] - freq_buf[major][0]) < 1.0f &&
                        fabsf(freq_buf[n][1] - freq_buf[major][1]) < 1.0f &&
                        fabsf(freq_buf[n][2] - freq_buf[major][2]) < 1.0f) vals[cnt++] = x_buf[n][j];
                x_mean[j] = TrimmedMean(vals, cnt);
            }
            // 频率取主组代表值(已规整)
            dw_f1 = freq_buf[major][0];
            dw_f2 = freq_buf[major][1];
            dw_f3 = freq_buf[major][2];
        }

        // dw_p1/2/3注释: LS拟合用fit_phase, 不依赖FFT相位
        // dw_p1 = fft8192_out_phase[0];
        // dw_p2 = fft8192_out_phase[1];
        // dw_p3 = fft8192_out_phase[2];

        k2 = (dw_f1 > 0.01f) ? (uint8_t)(dw_f2 / dw_f1) : 0;
        k3 = (dw_f1 > 0.01f) ? (uint8_t)(dw_f3 / dw_f1) : 0;
        sum_amp = MAX(dw_u1 + dw_u2 + dw_u3, 10.0f);

        switch (F_mode) {
        case 0: // 时域
            /*
                if (wave_src == 0) {
                    // 基础数字合成模式: 基波+谐波叠加(无相位)
                    DWIN_UpdateWaveBasic();
                }
                else if (wave_src == 1) {
                    // 直通模式: ADC数据直接映射
                    DWIN_UpdateWaveDirect();
                }
                else if (wave_src == 2) {
                    // 带相位合成模式: 使用FFT检测的频率+幅度+相位
                    DWIN_UpdateWaveSynth();
                }
                else {
                    // 最小二乘拟合模式: 时域ADC数据最优拟合
                    DWIN_UpdateWaveFit();
                }
            */

            DWIN_UpdateWaveFit();
            //发送波形
            DWIN_SendWave(0x04, wave_data, WAVE_POINTS);

            // LS拟合模式: 显示测量参数
            DWIN_SendText(TEXT1_ADDR, "Upp  = %3.1f mV", dw_upp);
            DWIN_SendText(TEXT2_ADDR, "Urms = %3.1f mV", dw_urms);
            DWIN_SendText(TEXT3_ADDR, "F1   = %3.1f kHz", dw_f1 / 1000.0f);

            break;
        case 1: // 频域
            if (MAX(k2, k3) <= 10) // 根据最大的谱线位置判别频谱图的美观的谱线间距
            {
                numx = 10;
                deltax = 70;
            }
            else if (MAX(k2, k3) <= 20) {
                numx = 20;
                deltax = 35;
            }
            else {
                numx = 50;
                deltax = 15;
            }
            /*
                        for (uint16_t i = 0; i < WAVE_POINTS; i++) {
                            if (i == deltax - 1) {
                                wave_data[i++] = 0;
                                wave_data[i++] = dw_u1 / sum_amp * WAVE_AMPL;
                                wave_data[i] = 0;
                            }
                            else if (i == k2 * deltax - 1) {
                                wave_data[i++] = 0;
                                wave_data[i++] = dw_u2 / sum_amp * WAVE_AMPL;
                                wave_data[i] = 0;
                            }
                            else if (i == k3 * deltax - 1) {
                                wave_data[i++] = 0;
                                wave_data[i++] = dw_u3 / sum_amp * WAVE_AMPL;
                                wave_data[i] = 0;
                            }
                            else wave_data[i] = 20;
                        }
            */
            for (uint16_t i = 0; i < WAVE_POINTS; i++) {
                if (i % deltax == 0) {
                    if (i == deltax) {
                        wave_data[i] = 2 * dw_u1 / sum_amp * WAVE_AMPL;
                    }
                    else if (i == k2 * deltax) {
                        wave_data[i] = 2 * dw_u2 / sum_amp * WAVE_AMPL;
                    }
                    else if (i == k3 * deltax) {
                        wave_data[i] = 2 * dw_u3 / sum_amp * WAVE_AMPL;
                    }
                    else {
                        wave_data[i] = 40;
                    }
                }
                else {
                    wave_data[i] = 20;
                }

            }
            DWIN_SendWave(0x04, wave_data, WAVE_POINTS);

            // 发送频域数据
            DWIN_SendText(TEXT1_ADDR, "%3.1f kHz , %3.1f mV", dw_f1 / 1000.0f, dw_u1);
            DWIN_SendText(TEXT2_ADDR, "%3.1f kHz , %3.1f mV", dw_f2 / 1000.0f, dw_u2);
            DWIN_SendText(TEXT3_ADDR, "%3.1f kHz , %3.1f mV", dw_f3 / 1000.0f, dw_u3);

            break;
        }
    }
}
