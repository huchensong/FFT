#include<fstream>//读文件用于测试
#include"../src/fft.h"
using namespace std;
//生成w table
int main(){
    std::ofstream outfile("/home/g/cpp_test/script/w_table_4096.txt", std::ios::out);
    outfile<<"Complex base_8_w[4096]={"<<endl;
    int num_coloum    = 16;
    int point_num       = 4096;
    int num_row = point_num/num_coloum; 
    for(int k = 0;k<num_row;k++){
    for(int i = 0; i<num_coloum;i++){
        outfile<<"Complex"<<calculate_w(1,i+16*k,4096)<<",";
    }
    outfile<<"\\"<<endl;
    }
    outfile<<"};"<<endl;
}