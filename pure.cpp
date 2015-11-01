#include "pure.h"

#include <iostream>
#include <iomanip>
using namespace std;




//Ѱ����ת�У�
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

//�ж���ת�����Ƿ�����ֵ������У����ص�һ���������±꣬���û�У�����-1��
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
//Ѱ����תԪ����������ת�У�
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

//����������⣻
//iCon, iVar, Origin  ���Լ�����̸������ܱ�����������������
void DZ(double pure[][MAXCol], int iCon, int iVar, int Origin)
{
	int i, j;
	int flag = -1;
	int row, col;
	float temp;
	float tp;

	while (1)
	{
		//�жϼ��ϵ���ǲ��Ƿ�����
		for (j = 0; j < iVar; j++)
			if (pure[iCon][j] > 0)
				flag = 1;
		if (flag == -1)
		{
			for (j = 0; j < iVar; j++)
				if (pure[iCon][j] < 0)
					flag++;
			if (flag == (Origin - 1))//������ϵ��Ϊ��ֵ�ĸ�������ԭʼ�����ĸ���ʱ��ֵΨһ��
			{
				printf("LP�����Ψһ����СֵΪ��%f\n", pure[iCon][iVar]);
				break;
			}
			else
			{
				printf("LP����ⲻΨһ����СֵΪ:%f\n", pure[iCon][iVar]);
				break;
			}
		}
		else
		{
			//Ѱ����ת��
			col = maxCol(pure, iCon, iVar);
			//�����ת�����г��˼��ϵ����������Ϊ�������������޽⣬����ѭ����
			if (isUpZero(pure, iCon, iVar, col) == -1)
			{
				printf("LP�����޽�!\n");
				break;
			}
			else
			{
				row = minRow(pure, iCon, iVar, col);
				//��λ����ת��
				printf("��תԪΪ��[%d,%d] %f\n", row, col, pure[row][col]);
				pure[row][iVar+1] = col;
				
				if (pure[row][col] != 1)
				{
					tp = pure[row][col];
					for (i = 0; i <= iVar + 1; i++)
						pure[row][i] /= tp;
				}
				//���ȱ任
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
				cout << "B�任��ĵ����α�Ϊ��" << endl;
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


//����ǵ�������  ʹ�ö��������
// iConԼ�����̸���   ��Ϊ�������У�����Լ�����̵ĸ���
// iVar�ܱ�������
// Originԭʼ��������
void UDZ(double pure[][MAXCol], int iCon, int iVar, int Origin)
{
	int i, j;
	int row, col;
	int flag = -1;
	float temp;
	float tp;
	//���K��ֵ������ô������
    //���ں���������������Լ�����̲����㵽��λ�����У����������Ҫ��2
	int k = iVar - iCon;
	while (1)
	{
		//�жϼ��ϵ���ǲ��Ƿ�����
		for (j = 0; j < iVar; j++)
			if (pure[iCon + 1][j] > 0)
				flag = 1;
		if (flag == -1)
		{
			if (pure[iCon + 1][iVar] == 0) //������ϵ��ȫΪ����ʱ������СֵΪ0��������תΪ�������⣬��������任��
			{
				for (i = 0; i <= iCon; i++)
				{
					pure[i][k]   = pure[i][iVar];
					pure[i][k+1] = pure[i][iVar+1];
				}

				//��������㷨ǰ������Ҫ��һ�����飬������Ҫ��֮ǰ���˹��������޳���
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
				//�޳��ɾ���
				
				DZ(pure, iCon, k, Origin);
				break;
			}
			else
			{
				printf("11LP�����޽�!\n");
				break;
			}
		}
		else
		{
			//Ѱ����ת��;
			col = maxCol(pure, iCon + 1, iVar);
			//�����ת�����г��˼��ϵ����������Ϊ�������������޽⣬����ѭ����
			if (isUpZero(pure, iCon + 1, iVar, col) == -1)
			{
				printf("LP�����޽�!\n");
				break;
			}
			row = minRow(pure, iCon, iVar, col);
			//��λ����ת��
			printf("��תԪΪ��[%d,%d] %f\n", row, col, pure[row][col]);
			pure[row][iVar+1] = col;

			if (pure[row][col] != 1)
			{
				tp = pure[row][col];
				for (i = 0; i < iVar + 1; i++)
					pure[row][i] /= tp;
			}
			//���ȱ任
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
			cout << "A�任��ĵ����б�Ϊ��" << endl;
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


