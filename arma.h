#ifndef __ARMA_H
#define __ARMA_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double MATRIX_TYPE;

typedef struct _Matrix
{
	/*Store Matrix
	存储矩阵*/
	int row;//行
	int column;//列
	MATRIX_TYPE *data;//矩阵指针
} Matirx;

float GetArmaCmathPredict(float *arinbuf,float *mainbuf,int p,int q,int sum);

void mattr(double *a,double *b,int m,int n);//通用矩阵转置
int minv(double *a,int n);//实方矩阵求逆
void rmmult(double *rm,double *a,double *b,int n,int m,int l);//兼容矩阵相乘

#endif
