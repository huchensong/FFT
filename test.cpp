#include <iostream>
#include <stdio.h>
#include <bits/stdc++.h>
using namespace std;
int n,a,b; 
int prime(int x)//这个函数是判断一个数是不是质数的
{	
	if(x==1)return 0;
	if(x==2)return 1;
	for(int i=2;i<=sqrt(x);i++)//i的初始值为什么是2以及为什么停止到sqrt（x）请看下面解析
	{
		if(x%i==0)	return 0;
	}
}



int main(){
	// int c=0;
	// while(cin>>n){
	// for(int i=n-4;i>=2;i--)
	// {
	// 	//printf("i=%d",i);
	// 	c=0;
	// 	if(prime(i)){
	// 		for(int k=n-i-1;k>1;){
	// 			//printf("k=%d",k);
	// 			k--;
	// 			if(prime(k)){
	// 				if(prime(n-k-i)){
	// 					c=3;
	// 					//printf("n-k-i=%d",n-k-i);
	// 				cout<<n-k-i<<" "<<k<<" "<<i<<endl;
	// 				break;
	// 				}
	// 			}
	// 		}
	// 		if(c==3)break;
	// 	}
	// }
	// }
	// return 0;
	cin>>n;
	for(int i=n-2;i>=n/2;i--){
		if(prime(i)){
			if(prime(n-i)){
				cout<<n<<"="<<i<<"+"<<n-i<<endl;
			}
		}
	}
}

