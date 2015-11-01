#ifndef _PURE_
#define _PURE_


#define MAXRow 300
#define MAXCol 300


/*单纯形法，典则问题求解
iCon: 约束方程个数
iVar: 总变量个数
Origin: 原始变量个数
*/
void DZ(double pure[][MAXCol], int iCon, int iVar, int Origin);
void UDZ(double pure[][MAXCol], int iCon, int iVar, int Origin);




#endif