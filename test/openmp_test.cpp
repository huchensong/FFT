 #include <iostream>
 #include <omp.h>   // openmp库
 #include<thread>
 #include<queue>
 #include<cmath>
 using namespace std;
 int main(){
 int sum = 0;
   auto start_time = std::chrono::high_resolution_clock::now();
   int cnt_ans;
   int vi[100];
   int len = 100;

     #pragma omp parallel for reduction(+:cnt_ans) default(shared)  num_threads(10)
        for(int i=0; i<len; i++)
    {
        auto &e = v_i[i]; 
        int t = 0; 

        if ( i < len / 2)
        {
            t = pow(e, 2); 
        }   
        else 
        {
            t = (int)sqrt(e); 
        }
        if ( t % 2 == 1)
        {
            cnt_ans += 1;
        }
    }
     auto end_time = std::chrono::high_resolution_clock::now();
     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
     cout << "Time taken by openmp: " << duration.count() << " microseconds" <<endl;
     cout << sum << endl;
//       start_time = std::chrono::high_resolution_clock::now();
//      for(int i=0; i<1000; i++)
//        {
//            sum +=  i;
//        }

//      end_time = std::chrono::high_resolution_clock::now();
//       duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
//      cout << "Time taken by no openmp: " << duration.count() << " microseconds" <<endl;
//  
    }

// g++ openmp_test.cpp -pthread -g -fopenmp -o openmp_test.out
