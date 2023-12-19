#include<thread>
#include<iostream>
#include<queue>

using namespace std;

void factorial(int &n){
    //计算n！
    cout<<n<<" fac"<<endl;
    int result=1;
    for(int i=1;i<=n;i++){
        result = result * i;
    }
    n = result;
    cout <<n<<endl;
}

void multiadd(int &n){
    //计算1+2+..+n
    cout<<n<<" add"<<endl;
    int a = 1;
    int result = 1;
    for(int i=1; i<=n; i++ ){
        result = result + i;
    }
    n = result;
    cout<<n<<endl;
}

void thread_test(){
    int a = 10;
    int b = 4000;
    auto start_time = std::chrono::high_resolution_clock::now();
    thread t1(factorial,std::ref(a));
    thread t2(multiadd,std::ref(b)) ;
    t1.join();
    t2.join();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    cout << "Time taken by function thread_test : " << duration.count() << " microseconds" <<endl;
}

void no_thread_test(){
    int a = 10;
    int b = 4000;
    auto start_time = std::chrono::high_resolution_clock::now();
    factorial(std::ref(a));
    multiadd(std::ref(b));
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    cout << "Time taken by function no_thread_test: " << duration.count() << " microseconds" <<endl;
}
int main(){
    thread_test();
    no_thread_test();
    return 0;
}