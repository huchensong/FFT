#include"fft.h"
#include<iostream>
#include<fstream>//读文件用于测试
#include<thread>
#include<iostream>
#include<queue>
#include<string>
using namespace std;
#define test_times 10
#define fft_length 512
#define prinr2file 0


int main(){
    //-------------读txt文件，准备测试数据-----------
    //如果test_times = 10 则读入10个txt文件的数据
    int16_t in_data[test_times][100000];
    int16_t out_result[test_times][100000];
    for(int i=0;i<test_times;i++){
        int data_length=0;
        ifstream  infile("/home/g/cpp_test/fft_data/fft_oai_in"+std::to_string(i)+".txt",ios::in);
        while(infile)infile>>in_data[i][data_length++];
    }
    //一半是实部一半是虚部
    // #define print_in_out
    #ifdef print_in_out  
    cout<<"len = "<<fft_length<<endl;
    cout<<"原始数据："<<endl;
    for(int i=0;i<2*fft_length;i++)cout<<in_data[i]<<' ';
    cout<<endl;
    #endif
    cout<<"---------------------------"<<endl;
    //--------------准备运行参数，进行FFT---------------------------
    int8_t base     = 8;
    int8_t in_type  = 0;
    int8_t out_type = 1;
    // if(base==2)out_type = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel for  num_threads(10)
    //-------------------运行FFT/dft------------------------------
    for(int test_id = 0; test_id<test_times; test_id++){
        fft_hcs(in_data[test_id],out_result[test_id],fft_length,base,in_type,out_type);
    }
    // for(int test_id = 0; test_id<test_times; test_id++){
    //     dft(in_data[test_id],out_result[test_id],fft_length,out_type);
    // }
    //-------------------end------------------------------

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    cout << "Time taken by function fft_hcs_base"+std::to_string(base)+"_version_1.0.1 with openmp: "<< duration.count()/test_times << " microseconds" <<endl;

    // dft_measure_time(in_data,out_result,fft_length,out_type);
    // fft_hcs_measure_time(in_data,out_result,fft_length,base,in_type,out_type);
    //-------------将结果写入文件-----------------------------------
    if(prinr2file){
    std::ofstream outfile("/home/g/cpp_test/fft_data/fft_oai_out_1_base"+std::to_string(base)+"_cpp.txt", std::ios::out);
    for (int i = 0; i < fft_length*2; i++) {
        outfile << out_result[i] << endl;
    }
    outfile.close(); 
    }
} 
//normal
// g++ fft.cpp -o fft.out
//with openmp 
// g++ fft.cpp -pthread -g -fopenmp -o fft.out