#include<thread>
#include<iostream>
#ifndef _FFT_H
#define _FFT_H

#include<queue>
#include <cmath>//指数，pi等
#include <complex>//支持复数运算
#include<fstream>//读文件用于测试
#include <omp.h>   // openmp库
using namespace std;
//头文件中不要放函数的定义，声明即可
//默认命名空间声明不要放在头文件，using namespace std;等应放在.cpp中，在 .h 文件中使用 std::string
typedef std:: complex<double> Complex;
Complex calculate_w(int16_t n,int16_t k,int16_t N)                                                                        ;
float   logM_n(int16_t M, int16_t n)                                                                                      ;
int     reverse(int n,int l)                                                                                              ;
void    dft(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type)                                           ;
void    dft_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type)                              ; 
void    fft_hcs (int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type,int8_t out_type)                         ;
void    fft_hcs_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type,int8_t out_type);
void    base2(int16_t i,int16_t data_gap,Complex w,int16_t* in_buf)                                                         ;
void    base8(int16_t i,int16_t data_gap,Complex* w,int16_t k,int16_t* in_buf)                                                         ;
void    stage_routine(int8_t stage,int16_t length,int8_t base,int16_t* in_buf,int16_t* out_buf)                           ;
int8_t  bitwidth(int16_t fft_length);

#endif
