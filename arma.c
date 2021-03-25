#include "arma.h"

static Matirx y = {0, 0, NULL}; //x y为计算所用矩阵
static Matirx x = {0, 0, NULL};

static Matirx *Matirx_creat(int row, int column);
static int M_free(Matirx *_mat);
static Matirx *CmathLeastSquares(void);
static void ArmaCmathInit(float *arbuf, float *mabuf, int p, int q, int sum);
static MATRIX_TYPE ArmaCmathPredict(Matirx *parameter, float *arinbuf, float *mainbuf, int p, int q, int sum);

/*
矩阵内存释放
*/
static int M_free(Matirx *_mat)
{
	free(_mat->data);
	free(_mat);
	return 0;
}
/*
使用ARMA模型预测下一个数据的值
*arinbuf ar模型部分输入，一般为实测时间序列
*mainbuf ma模型部分输入，一般为白噪声序列
p ar模型部分阶数，取值为大于等于1的整数
q ma模型部分阶数，取值为大于等于0的整数，当 q=0 时变为AR模型，预测依然有效
sum 单次预测所使用的样本数
*/
float GetArmaCmathPredict(float *arinbuf, float *mainbuf, int p, int q, int sum)
{
	Matirx *parameters = NULL;
	ArmaCmathInit(arinbuf, mainbuf, p, q, sum);
	parameters = CmathLeastSquares();
	/*
	printf("number of parameters:  %d\n", parameters->column);
	for (int i = 0; i < parameters->column; i++)
	{
		printf("parameters[%d] = %f\n", i, parameters->data[i]);
	}*/

	//预测下一个数据的值
	float predict_value = (float)ArmaCmathPredict(parameters, arinbuf, mainbuf, p, q, sum);
	//printf("predict value= %f\n", predict_value);
	M_free(parameters);
	free(x.data);
	free(y.data);
	return predict_value;
}
static MATRIX_TYPE ArmaCmathPredict(Matirx *parameter, float *arinbuf, float *mainbuf, int p, int q, int sum)
{
	MATRIX_TYPE s = 0;
	for (int i = 0; i < p; i++)
	{
		s += parameter->data[i] * arinbuf[sum - i - 1];
	}
	for (int i = 0; i < q; i++)
	{
		s += parameter->data[i + p] * mainbuf[sum - i - 1];
	}
	return s;
}
/*
ARMA模型计算初始化
*/
static void ArmaCmathInit(float *arbuf, float *mabuf, int p, int q, int sum)
{
	y.row = 1;
	y.column = sum - p - q;
	y.data = (MATRIX_TYPE *)malloc((y.row * y.column) * sizeof(MATRIX_TYPE));

	x.row = sum - p - q;
	x.column = p + q;
	x.data = (MATRIX_TYPE *)malloc((x.row * x.column) * sizeof(MATRIX_TYPE));
	/*printf("x.numRows=%d x.numCols%d\n", x.row, x.column);
		printf("\n");
		for(int j=0;j<sum;j++)
		{
			printf(" %f\t",arbuf[j]);
		}
		printf("\n\n");

		for(int k=0;k<sum;k++)
		{
			printf(" %f\t",mabuf[k]);
			
		}
		
		printf("\n\n");*/
	for (int i = 0; i < x.row; i++)
	{
		for (int j = 0; j < p; j++)
		{
			x.data[(p + q) * i + j] = arbuf[p + q - j + i - 1]; //printf("%d %f\t",(p+q)*i+j,x.data[(p+q)*i+j]);
		}

		for (int k = 0; k < q; k++)
		{
			x.data[(p + q) * i + p + k] = mabuf[p + q - k + i - 1]; //printf("%d %f\t",(p+q)*i+p+k,x.data[(p+q)*i+p+k]);
		}
		y.data[sum - p - q - 1 - i] = arbuf[sum - i - 1];

		//printf("\n");
	}
	/*for(int i=0;i<x.row;i++)
		{
			printf(" %f\n",y.data[i]);
		}
*/
}
/*
矩阵创建
*/
static Matirx *Matirx_creat(int row, int column)
{
	Matirx *p = malloc(sizeof(Matirx));
	p->row = row;
	p->column = column;
	p->data = malloc(sizeof(MATRIX_TYPE) * row * column);
	return p;
}
/**
 *最小二乘法求参数
 */
static Matirx *CmathLeastSquares(void)
{
	Matirx *y_t = Matirx_creat(y.column, y.row);
	mattr(y_t->data, y.data, y.row, y.column);
	Matirx *x_t = Matirx_creat(x.column, x.row);
	mattr(x_t->data, x.data, x.row, x.column);
	Matirx *m_x_t = Matirx_creat(x_t->row, x.column);
	rmmult(m_x_t->data, x_t->data, x.data, x_t->row, x_t->column, x.column);
	minv(m_x_t->data, m_x_t->row);
	Matirx *parameter1 = Matirx_creat(m_x_t->row, x_t->column);
	rmmult(parameter1->data, m_x_t->data, x_t->data, m_x_t->row, m_x_t->column, x_t->column);
	M_free(m_x_t);
	M_free(x_t);
	Matirx *parameter2 = Matirx_creat(parameter1->row, y_t->column);
	rmmult(parameter2->data, parameter1->data, y_t->data, parameter1->row, parameter1->column, y_t->column);
	M_free(parameter1);
	M_free(y_t);
	Matirx *parameter3 = Matirx_creat(parameter2->column, parameter2->row);
	mattr(parameter3->data, parameter2->data, parameter2->row, parameter2->column);
	M_free(parameter2);
	return parameter3;
}
