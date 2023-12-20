#include"fft.h"
#include<iostream>
#include<fstream>//读文件用于测试
#include<thread>
#include<iostream>
#include<queue>
#include<string>
using namespace std;
int main(){
    //-------------读txt文件，准备测试数据-----------
    ifstream  infile("/home/g/cpp_test/fft_data/fft_oai_in1.txt",ios::in);
    int16_t in_data[100000];
    int16_t out_result[100000];
    int data_length=0;
    while(infile)infile>>in_data[data_length++];
    data_length-=1;
    //一半是实部一半是虚部
    int fft_length = 512;
    data_length    = fft_length*2;
    // #define print_in_out
    #ifdef print_in_out  
    cout<<"len = "<<fft_length<<endl;
    cout<<"原始数据："<<endl;
    for(int i=0;i<2*fft_length;i++)cout<<in_data[i]<<' ';
    cout<<endl;
    #endif
    cout<<endl<<"---------------------------"<<endl;
    //--------------准备运行参数，进行FFT---------------------------
    int8_t base     = 2;
    std::ofstream outfile("/home/g/cpp_test/fft_data/fft_oai_out_1_base"+std::to_string(base)+"_cpp.txt", std::ios::out);
    int8_t in_type  = 0;
    int8_t out_type = 1;
    // if(base==2)out_type = 0;
    int test_times = 10;

    auto start_time = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel for  num_threads(10)
    //-------------------运行FFT/dft------------------------------
    for(int test_id = 0; test_id<test_times; test_id++){
        fft_hcs(in_data,out_result,fft_length,base,in_type,out_type);
    }
    // for(int test_id = 0; test_id<test_times; test_id++){
    //     dft(in_data,out_result,fft_length,out_type);
    // }
    //-------------------end------------------------------

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    cout << "Time taken by function fft_hcs_base8_version_1.0.0 without openmp: "<< duration.count()/test_times << " microseconds" <<endl;

    // dft_measure_time(in_data,out_result,fft_length,out_type);
    // fft_hcs_measure_time(in_data,out_result,fft_length,base,in_type,out_type);
    //-------------将结果写入文件-----------------------------------
    for (int i = 0; i < data_length; i++) {
        outfile << out_result[i] << endl;
    }
    outfile.close(); 
} 
//normal
// g++ fft.cpp -o fft.out
//with openmp 
// g++ fft.cpp -pthread -g -fopenmp -o fft.out