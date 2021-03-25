#include "arma.h"
#include "update_buf.h"

/*
使用ARMA或AR模型预测数据，本测试通过最小均方误差原则找到模型的最优参数，包括p，q，n
p为AR模型部分阶数，取值为大于等于1的整数
q为MA模型部分阶数，取值为大于等于0的整数，当 q=0 时变为AR模型，预测依然有效
n为单次预测所使用的样本数
*/
//所分析数据的最大长度
#define data_sum 400

double fangcha_min;//方差的默认最小值
int n_best, p_best, q_best;
int main()
{
    freopen("out.xls", "w", stdout);
    for (int n = 20; n <= 200; n += 20)//单次预测使用的样本数变化
    {
        predict_temp.buf_temp = (float *)malloc(n * sizeof(float));
        if (predict_temp.buf_temp == NULL)
            return 0;
        predict_temp.e_data = (float *)malloc(n * sizeof(float));
        if (predict_temp.e_data == NULL)
            return 0;
        fangcha_min = 100.;//方差的默认最小值
        for (int p = 1; p <= 5; p++)//AR模型部分的阶数变化
            for (int q = 0; q <= 5; q++)//MA模型部分的阶数变化
            {
                double fangcha = 0.;
                for (int i = 0; i < n; i++)
                {
                    update_bufs(&predict_temp, sensor_down[i], sensor_err[i], n);//初始化计算过程的bufs
                    //printf("%f %f \n",  predict_temp.buf_temp[n-1],  predict_temp.e_data[n-1]);
                }
                for (int i = n; i < data_sum; i++)
                {
                    predict_temp.predict_position = GetArmaCmathPredict(predict_temp.buf_temp, predict_temp.e_data, p, q, n);
                    predict_temp.predict_err_arma = sensor_down[i] - predict_temp.predict_position;
                    update_bufs(&predict_temp, sensor_down[i], predict_temp.predict_err_arma, n);
                    fangcha += predict_temp.predict_err_arma * predict_temp.predict_err_arma;
                }
                fangcha /= (float)(data_sum - n);
                //printf("%f\t",  predict_temp.predict_err_arma);
                //printf("n = %d\tp = %d\tq = %d\t%f\n", n, p, q, fangcha);
                if (fangcha < fangcha_min)//找到最优参数
                {
                    fangcha_min = fangcha;
                    n_best = n;
                    p_best = p;
                    q_best = q;
                }
            }
        printf("best result n = %d\tp = %d\tq = %d\tfangcha_min = %f\n", n_best, p_best, q_best, fangcha_min);
        // printf("%d\t%d\t%d\t%f\n", n_best, p_best, q_best, fangcha_min);
        free(predict_temp.buf_temp);
        free(predict_temp.e_data);
    }
    fclose(stdout);
    return 0;
}