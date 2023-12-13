#include <stdio.h>

int main() {
	int n ;
	scanf("%d", &n);
	printf("%d=", n);
	for (int i = 2; i <= n; i++) {
		while(n != i) { // 执行的条件必须是n与i不等，若相等则分解结束了
			if(n % i == 0) { // 若能整除则i为n的因子之一
				printf("%d*", i); // 输出因子
				n = n / i; // 找到了一个因子i，则n/i缩小n继续寻找
			} else {
				break; // 不能整除则跳出本次循环，递增i进行下一轮
            }
		}
	}
	printf("%d\n", n); // 最后剩下的n不能整除i，所以它也为因子之一，所以最后输出
	return 0;
}
