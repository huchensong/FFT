#include<thread>
#include<iostream>
#include<queue>
#include <cmath>//指数，pi等
#include <complex>//支持复数运算
#include<fstream>//读文件用于测试
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

typedef std:: complex<double> Complex;

Complex calculate_w(int16_t n,int16_t k,int16_t N)                                               ;
float   logM_n(int16_t M, int16_t n)                                                             ;
int     reverse(int n,int l)                                                                     ;
void    fft_hcs (int16_t* in_buf,int16_t* out_buf,int16_t fft_length, int8_t base,int8_t in_type);
void    base2(int16_t i,int16_t i_resp,Complex w,int16_t* in_buf)                                ;
void    stage_routine(int8_t stage,int16_t length,int8_t base,int16_t* in_buf,int16_t* out_buf)  ;
int8_t  bitwidth(int16_t fft_length);
//第一种计算旋转因子的方式，直接使用定义计算,计算量较大,计算Wn_k = exp(-2*pi/N*n*k) 其中pi = acos(-1)
Complex calculate_w(int16_t n,int16_t k,int16_t N){
    //计算Wn_k = exp(-2*pi/N*n*k) 其中pi = acos(-1)
    Complex c;
    c.real(cos(2*acos(-1)*n*k/N));
    c.imag(sin(-2*acos(-1)*n*k/N));
    return c;
}

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
int8_t  bitwidth(int16_t fft_length){
    int8_t bit_width = 0;
    while(fft_length){
        fft_length = fft_length>>1;
        bit_width++;
    }
    return bit_width;
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
            cout<<"outbuf"<<reverse(2*i,bitwidth(fft_length)-1)<<"=inbuf"<<2*i<<endl;
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

void base2(int16_t i,int16_t i_resp,Complex w,int16_t* in_buf){
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
}


void stage_routine(int8_t stage,int16_t fft_length,int8_t base,int16_t* in_buf,int16_t* out_buf){
    //stage 决定取数据的间隔 ,要求length是2的幂次,
    //if length=2048 stage = 0 ,base =2 ,then  取数据间隔data_gap = length/base/(2^stage) = 1024 此时做的是模pow(base,stage+1)的FFT
    //疑问：2的幂次，使用pow函数与直接移位会有时延差别吗
    //解决疑问方法1：查看反汇编代码
    //解决疑问方法2：用循环来实际测试验证
    //验证结果：
    int16_t data_gap = fft_length/base/pow(2,stage);
    cout<<" fft_length ="<<fft_length<<" stage ="<<(int)stage<<" datagap = "<<data_gap<<endl;
    //分block进行运算，先计算存在几个block，每个block的起始地址相隔的距离
    int16_t block_num    = pow(base,stage)     ;
    int16_t block_length = fft_length/block_num;
    int16_t block_cnt    = block_length/base   ;
    int16_t cnt = 0;
    for(int k=0; k<block_num;k++){
        for(int cnt =0;cnt<block_cnt;cnt++){
            int N = pow(base,stage+1);
            int16_t base_addr = block_length*k;
            switch (base)
            {
                case 2:
                    base2(base_addr + cnt,base_addr + cnt + data_gap,calculate_w(1,reverse(k,stage),N),in_buf);
                    //找不到旋转因子计算的规律....这一版就写到这,这一句有问题
                    break;
                default:
                    break;
            }
        }
    } 
}


int main(){
    //-------------读txt文件，准备测试数据-----------
    ifstream  infile("fft_data/fft_oai_in1.txt",ios::in);
    int16_t in_data[100000];
    int16_t out_result[100000];
    int data_length=0;
    while(infile)infile>>in_data[data_length++];
    data_length-=1;
    //一半是实部一半是虚部
    int fft_length = data_length/2;
    #ifdef print_in_out  
    cout<<"len = "<<fft_length<<endl;
    cout<<"原始数据："<<endl;
    for(int i=0;i<fft_length;i++)cout<<in_data[i]<<' ';
    #endif
    cout<<logM_n(8,4096);
    cout<<endl<<"---------------------------"<<endl;
    //--------------准备运行参数，进行FFT---------------------------
    int8_t base     = 2;
    int8_t in_type  = 0;
    int8_t out_type = 0;
    fft_hcs(in_data,out_result,fft_length,base,in_type,out_type);
    #ifdef print_in_out  
    for(int i=0;i<fft_length;i++)cout<<in_data[i]<<' ';
    #endif
    //-------------将结果写入文件-----------------------------------
    std::ofstream outfile("fft_data/fft_oai_out_1_cpp.txt", std::ios::out);
    for (int i = 0; i < data_length; i++) {
        outfile << out_result[i] << endl;
    }
    cout<<(int)bitwidth(2)<<endl;
    cout<<(int)bitwidth(2048)<<endl;
    cout<<(int)bitwidth(1024)<<endl;
    // outfile.close(); 
} 