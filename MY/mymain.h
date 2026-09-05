//zzy：此c及h可废弃不用，内容已经整合至DWIN
#ifndef MY_MAIN_H
#define MY_MAIN_H

#include "adcdac.h"


// 工作模式
typedef enum {
    MODE_IDLE = 0,    //空闲模式
    MODE_PASS = 1, //直通模式
    MODE_SCOP = 2,     //示波模式
} WorkMode;

// 全局变量
extern WorkMode mode;

// 主循环处理
void mymain_process(void);


#endif /* MY_MAIN_H */
