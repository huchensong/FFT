
#include <thread>
#include <future>
#include <numeric>
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
 
double accumulate(int min, int max)
{
	double sum = 0;
	for (int i = min; i <= max; ++i)
	{
		sum += sqrt(i);
	}
    std::cout<<"here"<<std::endl;
	return sum;
}
 
 
double concurrent_task(int min, int max)
{
	//每个线程的执行结果存放在容器中
	std::vector<std::future<double>> results;
	unsigned concurrent_count = std::thread::hardware_concurrency();
	min = 0;
	for (int i = 0; i < concurrent_count; i++)
	{
		std::packaged_task<double(int, int)> task(accumulate); //产生一个未就绪的共享状态
		results.push_back(task.get_future());
 
		int range = max / concurrent_count * (i + 1); //任务平均分配到各个线程中
        std::cout<< concurrent_count<<std::endl;
		std::thread t(std::move(task), min, range); //通过新线程执行任务
		t.detach();
 
		min = range + 1;
	}
 
	std::cout << "threads create finish" << std::endl;
	double sum = 0;
	for (auto& r : results) {
		sum += r.get(); // 通过future获取每个任务的结果，即获取共享状态
	}
	return sum;
}
 
int main()
{
	auto start_time = std::chrono::steady_clock::now();
 
	double r = concurrent_task(1, 10e8);
 
	auto end_time = std::chrono::steady_clock::now();
 
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
 
	std::cout << "Concurrent task finish, " << ms << " ms consumed, Result: " << r << std::endl;

    start_time = std::chrono::steady_clock::now();
 
	r = accumulate(1, 10e8);
 
	end_time = std::chrono::steady_clock::now();
 
	ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
 
	std::cout << "Common task finish, " << ms << " ms consumed, Result: " << r << std::endl;

	return 0;

}