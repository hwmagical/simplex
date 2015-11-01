#include "pure.h"

#include <iostream>
#include <iomanip>
using namespace std;




//寻找旋转列；
int maxCol(double array[][MAXCol], int iCon, int iVal)
{
	int j;
	float max = array[iCon][0];
	int jmax = 0;
	for (j = 1; j < iVal; j++)
	{
		if (array[iCon][j] > max)
		{
			max = array[iCon][j];
			jmax = j;
		}
	}
	if (max > 0)
		return jmax;
	return -1;
}

//判断旋转列中是否有正值，如果有，返回第一个正数的下标，如果没有，返回-1；
int isUpZero(double array[][MAXCol], int iCon, int iVal, int col)
{
	int flag = -1;
	int i;
	for (i = 0; i < iCon; i++)
	{
		if (array[i][col] > 0)
		{
			flag = i;
			break;
		}
	}
	return flag;
}
//寻找旋转元，并返回旋转行；
int minRow(double array[][MAXCol], int iCon, int iVal, int col)
{
	int i;
	float min;
	int upZero = isUpZero(array, iCon, iVal, col);
	int flag = upZero;
	if (upZero != -1)
	{
		min = array[upZero][iVal] / array[upZero][col];
		for (i = upZero + 1; i < iCon; i++)
			if (array[i][col] > 0 && array[i][iVal] / array[i][col] < min)
			{
			flag = i;
			min = array[i][iVal] / array[i][col];
			}

	}
	return flag;

}

//解决典则问题；
//iCon, iVar, Origin  存放约束方程个数，总变量个数，变量个数
void DZ(double pure[][MAXCol], int iCon, int iVar, int Origin)
{
	int i, j;
	int flag = -1;
	int row, col;
	float temp;
	float tp;

	while (1)
	{
		//判断检测系数是不是非正；
		for (j = 0; j < iVar; j++)
			if (pure[iCon][j] > 0)
				flag = 1;
		if (flag == -1)
		{
			for (j = 0; j < iVar; j++)
				if (pure[iCon][j] < 0)
					flag++;
			if (flag == (Origin - 1))//当检验系数为负值的个数等于原始变量的个数时，值唯一；
			{
				printf("LP问题解唯一，最小值为：%f\n", pure[iCon][iVar]);
				break;
			}
			else
			{
				printf("LP问题解不唯一，最小值为:%f\n", pure[iCon][iVar]);
				break;
			}
		}
		else
		{
			//寻找旋转列
			col = maxCol(pure, iCon, iVar);
			//如果旋转列所有除了检测系数，其他均为非正，则问题无解，结束循环；
			if (isUpZero(pure, iCon, iVar, col) == -1)
			{
				printf("LP问题无解!\n");
				break;
			}
			else
			{
				row = minRow(pure, iCon, iVar, col);
				//单位化旋转行
				printf("旋转元为：[%d,%d] %f\n", row, col, pure[row][col]);
				pure[row][iVar+1] = col;
				
				if (pure[row][col] != 1)
				{
					tp = pure[row][col];
					for (i = 0; i <= iVar + 1; i++)
						pure[row][i] /= tp;
				}
				//初等变换
				for (i = 0; i <= iCon;)
				{
					if (i == row)
						i++;
					else
					{
						temp = pure[i][col];
						for (j = 0; j <= iVar; j++)
							pure[i][j] -= temp*pure[row][j];
						i++;
					}
				}
				
				/*
				cout << "B变换后的单纯形表为：" << endl;
				for (j = 0; j < iVar + 2; j++)
					cout << left << setw(8) << j;
				cout << endl;
				for (i = 0; i <= iCon; i++)
				{
					for (j = 0; j < iVar + 2; j++)
						cout <<left<<setw(8) << pure[i][j];
					cout << endl;
				}
				*/
				flag = -1;
			}
		}
	}
}


//解决非典则问题  使用二步法求解
// iCon约束方程个数   仅为二步法中，正真约束方程的个数
// iVar总变量个数
// Origin原始变量个数
void UDZ(double pure[][MAXCol], int iCon, int iVar, int Origin)
{
	int i, j;
	int row, col;
	int flag = -1;
	float temp;
	float tp;
	//这个K数值不能这么简单来算
    //由于航班问题的最后两个约束方程不计算到单位矩阵中，因此这里需要加2
	int k = iVar - iCon;
	while (1)
	{
		//判断检测系数是不是非正；
		for (j = 0; j < iVar; j++)
			if (pure[iCon + 1][j] > 0)
				flag = 1;
		if (flag == -1)
		{
			if (pure[iCon + 1][iVar] == 0) //当检验系数全为负数时，且最小值为0，则问题转为典则问题，以下做表变换；
			{
				for (i = 0; i <= iCon; i++)
				{
					pure[i][k]   = pure[i][iVar];
					pure[i][k+1] = pure[i][iVar+1];
				}

				//进入典则算法前，还需要做一件事情，就是需要把之前的人工变量都剔除掉
				for (i = 0; i < iCon; i++)
				{
					if (pure[i][k+1] == 0xffffffff)
					{
						for (j = 0; j < k; j++)
						{
							pure[i][j] = 0;
						}
					} 
				}
				for (j = 0; j < k; j++)
				{
					if (pure[iCon + 1][j] != 0)
					{
						for (i = 0; i < iCon + 1; i++)
						{
							pure[i][j] = 0;
						}
					}
				}
				//剔除干净了
				
				DZ(pure, iCon, k, Origin);
				break;
			}
			else
			{
				printf("11LP问题无解!\n");
				break;
			}
		}
		else
		{
			//寻找旋转列;
			col = maxCol(pure, iCon + 1, iVar);
			//如果旋转列所有除了检测系数，其他均为非正，则问题无解，结束循环；
			if (isUpZero(pure, iCon + 1, iVar, col) == -1)
			{
				printf("LP问题无解!\n");
				break;
			}
			row = minRow(pure, iCon, iVar, col);
			//单位化旋转行
			printf("旋转元为：[%d,%d] %f\n", row, col, pure[row][col]);
			pure[row][iVar+1] = col;

			if (pure[row][col] != 1)
			{
				tp = pure[row][col];
				for (i = 0; i < iVar + 1; i++)
					pure[row][i] /= tp;
			}
			//初等变换
			for (i = 0; i <= iCon + 1;)
			{
				if (i == row)
					i++;
				else
				{
					temp = pure[i][col];
					for (j = 0; j <= iVar; j++)
						pure[i][j] -= temp*pure[row][j];
					i++;
				}
			}
			
			/*
			cout << "A变换后的单纯行表为：" << endl;
			for (j = 0; j < iVar + 2; j++)
				cout << left << setw(8) << j;
			cout << endl;
			for (i = 0; i <= iCon + 1; i++)
			{
				for (j = 0; j < iVar + 2; j++)
					cout << left << setw(8) << pure[i][j];
				cout << endl;
			}
			*/
			flag = -1;
		}
	}
}

//hahahahaha


