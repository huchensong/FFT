#include<thread>
#include<iostream>
#include<queue>
#include <cmath>//指数，pi等
#include <complex>//支持复数运算
#include<fstream>//读文件用于测试
#include <omp.h>   // openmp库
#include <vector>
#include <future>
#include <chrono>
#include"fft.h"
#define OPENMP_ENABLE
#define pi 3.1415926
using namespace std;



/*
先不考虑多线程，实现FFT之后再考虑多线程分解任务
--------------设计思路----------------------
目标：给定数据，自动取相适应的2的幂次方长度进行FFT
version 1.0.0: void fft(int16_t* in_buf,int16_t* out_buf);

parameter in_type:
in_type==0 :in_buf1 以自然位序输入
in_type==1 :in_buf1 以自然倒位序输入

parameter out_type:
out_type==0 :out_buf1 以自然位序输出
out_type==1 :out_buf1 以自然倒位序输出

parameter base:
base==2 :基2
base==4 :基4
base==8 :基8

parameter fft_length: 决定fft长度，暂时只支持2的幂次


问题1：
旋转因子采取计算还是查表的方式？
按理说查表的方式精度应该够了

实现思路：

*/

// typedef std:: complex<double> Complex;

// Complex calculate_w(int16_t n,int16_t k,int16_t N)                                                                        ;
// float   logM_n(int16_t M, int16_t n)                                                                                      ;
// int     reverse(int n,int l)                                                                                              ;
// void    dft(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type)                                           ;
// void    dft_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type)                              ; 
// void    fft_hcs (int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type)                         ;
// void    fft_hcs_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type,int8_t out_type);
// void    base2(int16_t i,int16_t data_gap,Complex w,int16_t* in_buf)                                                         ;
// void    base8(int16_t i,int16_t data_gap,Complex w,int16_t* in_buf)                                                         ;
// void    stage_routine(int8_t stage,int16_t length,int8_t base,int16_t* in_buf,int16_t* out_buf)                           ;
// int8_t  bitwidth(int16_t fft_length);

//第一种计算旋转因子的方式，直接使用定义计算,计算量较大,计算Wn_k = exp(-2*pi/N*n*k) 其中pi = acos(-1)


Complex calculate_w(int16_t n,int16_t k,int16_t N){
    //计算Wn_k = exp(-2*pi/N*n*k) 其中pi = acos(-1)
    Complex c;
    c.real(cos(2*acos(-1)*n*k/N));
    c.imag(sin(-2*acos(-1)*n*k/N));
    return c;
}

Complex calculate_int_w(int16_t n,int16_t k,int16_t N){
    //计算Wn_k = exp(-2*pi/N*n*k) 其中pi = acos(-1)
    Complex c;
    c.real(cos(2*acos(-1)*n*k/N));
    c.imag(sin(-2*acos(-1)*n*k/N));
    return c;
}
//第二种计算旋转因子的方式，查表方式，memory花费较大，以空间换时间
//首先考虑二维数组，但由于旋转因子的稀疏性，会有很多没用到的空间，利用率太低
//再考虑类似字典的数据结构，hash表？
//有两种方式 一种是在开始之前init table，第二种是提前给初始值，用代码的方式初始化



float logM_n(int16_t M, int16_t n){
    return log2(n)/log2(M);
}

int reverse(int n,int l){
    int i=0,r=0,t;
    while(i<l){
       t=n&1;
       n=n>>1;
       r+=t*pow(2,l-1-i); 
       i++;
    }
    return r;
}
//返回fft_length的位长度 比如输入1024 返回 11
int8_t  bitwidth(int16_t fft_length){
    int8_t bit_width = 0;
    while(fft_length){
        fft_length = fft_length>>1;
        bit_width++;
    }
    return bit_width;
}

//根据定义式X[k]=sum{x(n)*Wnk}计算
void    dft(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type){
    for(int k=0;k<fft_length;k++){
        for(int n=0;n<fft_length;n++){
            out_buf[2*k]   += in_buf[2*n]*calculate_w(n,k,fft_length).real() - in_buf[2*n+1]*calculate_w(n,k,fft_length).imag(); 
            out_buf[2*k+1] += in_buf[2*n]*calculate_w(n,k,fft_length).imag() + in_buf[2*n+1]*calculate_w(n,k,fft_length).real();
        }
    }
}

void    dft_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length,int8_t in_type){
    auto start_time = std::chrono::high_resolution_clock::now();
    dft(in_buf,out_buf,fft_length,in_type);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Time taken by function dft: "<< duration.count() << " microseconds" <<std::endl;
}

void fft_hcs (int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type,int8_t out_type){
    //计算在此base下需要经过几个stage才能输出最终结果
    int8_t stage_num = logM_n(base,fft_length);
    for(int stage_id = 0; stage_id <stage_num; stage_id++){
        stage_routine(stage_id,fft_length,base,in_buf,out_buf);
    }
    // 将in_buf中的数据原封不动/调整位序后输入out_buf
    for(int i = 0; i<fft_length;i++){
        switch (out_type)
        {
        case 0:
            //out_type=0,自然位序输出
            // cout<<"outbuf"<<reverse(2*i,bitwidth(fft_length)-1)<<"=inbuf"<<2*i<<endl;
            out_buf[2*reverse(i,bitwidth(fft_length)-1)]   = in_buf[2*i]  ;
            out_buf[2*reverse(i,bitwidth(fft_length)-1)+1] = in_buf[2*i+1];
            break;
        case 1:
            out_buf[2*i]   = in_buf[2*i]  ;
            out_buf[2*i+1] = in_buf[2*i+1];
            break;
        default:
            out_buf[2*i]   = in_buf[2*i]  ;
            out_buf[2*i+1] = in_buf[2*i+1];
        }
    }
}
void fft_hcs_measure_time(int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type,int8_t out_type){
    auto start_time = std::chrono::high_resolution_clock::now();
    fft_hcs(in_buf,out_buf,fft_length,base,in_type,out_type);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    cout << "Time taken by function fft_base"<<base<<": "<< duration.count() << " microseconds" <<endl;
}
void base2(int16_t i,int16_t data_gap,Complex w,int16_t* in_buf){
    //这里使用data_gap是为了和后面的base8参数列表呼应，也可以使用data_gap，而不用传8个索引值
    int16_t i_resp = i + data_gap;
    //i对应的实部是in_buf[2*i], 虚部为 in_buf[2*i+1]
    int16_t re_temp        = in_buf[2*i]                                          ;
    int16_t im_temp        = in_buf[2*i+1]                                        ;
    int16_t b_product_w_re = in_buf[2*i_resp]*w.real()-in_buf[2*i_resp+1]*w.imag();
    int16_t b_product_w_im = in_buf[2*i_resp]*w.imag()+in_buf[2*i_resp+1]*w.real();
    //计算实部 a_re = a_re + b_re * w_re - b_im * w_im ;
    //计算虚部 a_im = a_im + b_im * w_re + b_re * w_im ;
    in_buf[2*i]        = re_temp + b_product_w_re;
    in_buf[2*i+1]      = im_temp + b_product_w_im;
    in_buf[2*i_resp]   = re_temp - b_product_w_re;
    in_buf[2*i_resp+1] = im_temp - b_product_w_im;
    //cout<<in_buf[0]<<" ";
}
// void base2(int16_t i,int16_t data_gap,Complex* w,int16_t* in_buf){
//     //这里使用data_gap是为了和后面的base8参数列表呼应，也可以使用data_gap，而不用传8个索引值
//     int16_t i_resp = i + data_gap;
//     //i对应的实部是in_buf[2*i], 虚部为 in_buf[2*i+1]
//     int16_t re_temp        = in_buf[2*i]                                          ;
//     int16_t im_temp        = in_buf[2*i+1]                                        ;
//     int16_t b_product_w_re = in_buf[2*i_resp]*w.real()-in_buf[2*i_resp+1]*w.imag();
//     int16_t b_product_w_im = in_buf[2*i_resp]*w.imag()+in_buf[2*i_resp+1]*w.real();
//     //计算实部 a_re = a_re + b_re * w_re - b_im * w_im ;
//     //计算虚部 a_im = a_im + b_im * w_re + b_re * w_im ;
//     in_buf[2*i]        = re_temp + b_product_w_re;
//     in_buf[2*i+1]      = im_temp + b_product_w_im;
//     in_buf[2*i_resp]   = re_temp - b_product_w_re;
//     in_buf[2*i_resp+1] = im_temp - b_product_w_im;
//     //cout<<in_buf[0]<<" ";
// }
// void base2(int16_t i,int16_t i_resp,Complex w,int16_t* in_buf){
//     //i对应的实部是in_buf[2*i], 虚部为 in_buf[2*i+1]
//     int16_t re_temp        = in_buf[2*i]                                          ;
//     int16_t im_temp        = in_buf[2*i+1]                                        ;
//     int16_t b_product_w_re = in_buf[2*i_resp]*w.real()-in_buf[2*i_resp+1]*w.imag();
//     int16_t b_product_w_im = in_buf[2*i_resp]*w.imag()+in_buf[2*i_resp+1]*w.real();
//     //计算实部 a_re = a_re + b_re * w_re - b_im * w_im ;
//     //计算虚部 a_im = a_im + b_im * w_re + b_re * w_im ;
//     in_buf[2*i]        = re_temp + b_product_w_re;
//     in_buf[2*i+1]      = im_temp + b_product_w_im;
//     in_buf[2*i_resp]   = re_temp - b_product_w_re;
//     in_buf[2*i_resp+1] = im_temp - b_product_w_im;
// }

void base8(int16_t i,int16_t data_gap,Complex* w,int16_t k,int16_t* in_buf){
    Complex minus_j(0,-1)           ;
    double sqrt_2   = sqrt(0.5)        ;
    Complex W8_1(sqrt_2,-sqrt_2) ;
    Complex W8_3(-sqrt_2,-sqrt_2);

    //--------------------分离实部虚部------------------

    Complex temp_0_0(in_buf[2*i]             ,in_buf[2*i+1]);
    Complex temp_1_0(in_buf[2*(i+data_gap)]  ,in_buf[2*((i+data_gap))+1]);
    Complex temp_2_0(in_buf[2*(i+2*data_gap)],in_buf[2*((i+2*data_gap))+1]);
    Complex temp_3_0(in_buf[2*(i+3*data_gap)],in_buf[2*((i+3*data_gap))+1]);
    Complex temp_4_0(in_buf[2*(i+4*data_gap)],in_buf[2*((i+4*data_gap))+1]);
    Complex temp_5_0(in_buf[2*(i+5*data_gap)],in_buf[2*((i+5*data_gap))+1]);
    Complex temp_6_0(in_buf[2*(i+6*data_gap)],in_buf[2*((i+6*data_gap))+1]);
    Complex temp_7_0(in_buf[2*(i+7*data_gap)],in_buf[2*((i+7*data_gap))+1]);
    // cout<<"temp_0_0="<<temp_0_0<<" ";
    // cout<<"temp_1_0="<<temp_1_0<<" ";
    // cout<<"temp_2_0="<<temp_2_0<<" ";
    // cout<<"temp_3_0="<<temp_3_0<<" ";
    // cout<<"temp_4_0="<<temp_4_0<<" ";
    // cout<<"temp_5_0="<<temp_5_0<<" ";
    // cout<<"temp_6_0="<<temp_6_0<<" ";
    // cout<<"temp_7_0="<<temp_7_0<<" ";
    // cout<<endl;
    //------------------与旋转因子相乘-----------------
    Complex temp_0_1        =temp_0_0*w[0];
    Complex temp_1_1        =temp_1_0*w[k];
    Complex temp_2_1        =temp_2_0*w[2*k];
    Complex temp_3_1        =temp_3_0*w[3*k];
    Complex temp_4_1        =temp_4_0*w[4*k];
    Complex temp_5_1        =temp_5_0*w[5*k];
    Complex temp_6_1        =temp_6_0*w[6*k];
    Complex temp_7_1        =temp_7_0*w[7*k];
    //-----------------stage 1 -----------------------
    temp_0_0  = temp_0_1 + temp_4_1;
    temp_4_0  = temp_0_1 - temp_4_1;
    temp_1_0  = temp_1_1 + temp_5_1;
    temp_5_0  = temp_1_1 - temp_5_1;
    temp_2_0  = temp_2_1 + temp_6_1;
    temp_6_0  = temp_2_1 - temp_6_1;
    temp_3_0  = temp_3_1 + temp_7_1;
    temp_7_0  = temp_3_1 - temp_7_1;
    //-----------------stage 2 -----------------------
    temp_0_1  = temp_0_0 + temp_2_0;
    temp_2_1  = temp_0_0 - temp_2_0;
    temp_1_1  = temp_1_0 + temp_3_0;
    temp_3_1  = temp_1_0 - temp_3_0;
    temp_4_1  = temp_4_0 + temp_6_0 * minus_j;
    temp_6_1  = temp_4_0 - temp_6_0 * minus_j;
    temp_5_1  = temp_5_0 + temp_7_0 * minus_j;
    temp_7_1  = temp_5_0 - temp_7_0 * minus_j;
    //--------------stage 3----------------------------
    temp_0_0  = temp_0_1 + temp_1_1;
    temp_1_0  = temp_0_1 - temp_1_1;
    temp_2_0  = temp_2_1 + temp_3_1 * minus_j;
    temp_3_0  = temp_2_1 - temp_3_1 * minus_j;
    temp_4_0  = temp_4_1 + temp_5_1*W8_1;
    temp_5_0  = temp_4_1 - temp_5_1*W8_1;
    temp_6_0  = temp_6_1 + temp_7_1*W8_3;
    temp_7_0  = temp_6_1 - temp_7_1*W8_3;
    //----------赋值给in_buf----------------------------
    in_buf[2*i]                = temp_0_0.real();
    in_buf[2*i+1]              = temp_0_0.imag();

    in_buf[2*(i+data_gap)]     = temp_4_0.real();
    in_buf[2*(i+data_gap)+1]   = temp_4_0.imag();

    in_buf[2*(i+2*data_gap)]   = temp_2_0.real();
    in_buf[2*(i+2*data_gap)+1] = temp_2_0.imag();

    in_buf[2*(i+3*data_gap)]   = temp_6_0.real();
    in_buf[2*(i+3*data_gap)+1] = temp_6_0.imag();

    in_buf[2*(i+4*data_gap)]   = temp_1_0.real();
    in_buf[2*(i+4*data_gap)+1] = temp_1_0.imag();

    in_buf[2*(i+5*data_gap)]   = temp_5_0.real();
    in_buf[2*(i+5*data_gap)+1] = temp_5_0.imag();

    in_buf[2*(i+6*data_gap)]   = temp_3_0.real();
    in_buf[2*(i+6*data_gap)+1] = temp_3_0.imag();

    in_buf[2*(i+7*data_gap)]   = temp_7_0.real();
    in_buf[2*(i+7*data_gap)+1] = temp_7_0.imag();
    // cout<<"temp_0_0="<<temp_0_0<<" ";
    // cout<<"temp_1_0="<<temp_1_0<<" ";
    // cout<<"temp_2_0="<<temp_2_0<<" ";
    // cout<<"temp_3_0="<<temp_3_0<<" ";
    // cout<<"temp_4_0="<<temp_4_0<<" ";
    // cout<<"temp_5_0="<<temp_5_0<<" ";
    // cout<<"temp_6_0="<<temp_6_0<<" ";
    // cout<<"temp_7_0="<<temp_7_0<<" ";
    // cout<<endl;
}

void stage_routine(int8_t stage,int16_t fft_length,int8_t base,int16_t* in_buf,int16_t* out_buf){
    //stage 决定取数据的间隔 ,要求length是2的幂次,
    //if length=2048 stage = 0 ,base =2 ,then  取数据间隔data_gap = length/base/(2^stage) = 1024 此时做的是模pow(base,stage+1)的FFT
    //疑问：2的幂次，使用pow函数与直接移位会有时延差别吗
    //解决疑问方法1：查看反汇编代码
    //解决疑问方法2：用循环来实际测试验证
    //验证结果：
    int16_t block_num    = pow(base,stage)     ;
    int16_t block_length = fft_length/block_num;
    int16_t data_gap = block_length/base;
    unsigned concurrent_count = std::thread::hardware_concurrency();
    //cout<<"using"<<concurrent_count<<"threads";
    // cout<<" fft_length ="<<fft_length<<" stage ="<<(int)stage<<" datagap = "<<data_gap<<endl;
    //分block进行运算，先计算存在几个block，每个block的起始地址相隔的距离
    int16_t block_cnt    = data_gap   ;
    for(int k=0; k<block_num;k++){
        for(int cnt =0;cnt<block_cnt;cnt++){
            int N = pow(base,stage+1);
            int16_t base_addr = block_length*k;
            // std::packaged_task<void(int16_t,int16_t,Complex*,int16_t,int16_t*)>task(base8);
            // std::thread t(std::move(task), base_addr + cnt,data_gap,base_8_w,64*(k&7)+(k&0b111111000),in_buf); //通过新线程执行任务
            // t.detach();
            // switch (base)
            // {
            //     case 2:
            //         base2(base_addr + cnt,data_gap,calculate_w(1,reverse(k,stage),N),in_buf);
            //         //  base2(base_addr + cnt,base_addr + cnt + data_gap,calculate_w(1,reverse(k,stage),N),in_buf);
            //         break;
            //     case 8:
            //         // cout<<"base_addr+cnt="<<base_addr + cnt<<endl;
            //         // cout<<"N=4096,n*k="<<64*(k&7)+(k&0b111111000)<<endl;
            //         // std::packaged_task<void(int16_t,int16_t,Complex*,int16_t,int16_t*)>task(base8);
            //         // std::thread t(std::move(task), base_addr + cnt,data_gap,base_8_w,64*(k&7)+(k&0b111111000),in_buf); //通过新线程执行任务
            //         // t.detach();
            //         base8(base_addr + cnt,data_gap,base_8_w,64*(k&7)+(k&0b111111000),in_buf);
            //         break;
            //     default:
            //         cout<<"base input unsupported!"<<endl;
            //         break;
            // }
        }
    } 
}
