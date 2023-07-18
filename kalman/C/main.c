#include "stdio.h"
static double Q = 0.18;  //过程噪声协方差
static double R = 0.54;    //观测噪声协方差（大）

//矩阵相乘
//rows1 表示 m1 的行数，cols1 表示 m1 的列数（也是 m2 的行数），cols2 表示 m2 的列数
//m1[rows1][cols1]
//m2[cols1][cols2]
//result[rows1][cols2]
void matrix3v3(double m1[][3], double m2[][3], double result[][3], int rows1, int cols1, int cols2)
{
    int i, j, k;

    for (i = 0; i < rows1; i++) {
        for (j = 0; j < cols2; j++) {
            result[i][j] = 0;
            for (k = 0; k < cols1; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}
//
void matrix3v1(double m1[][3], double m2[][1], double result[][1], int rows1, int cols1, int cols2)
{
    int i, j, k;

    for (i = 0; i < rows1; i++) {
        for (j = 0; j < cols2; j++) {
            result[i][j] = 0;
            for (k = 0; k < cols1; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

//转置矩阵3x3
void transm3x3(double matrix[][3], double result[][3], int rows, int columns)
{
    // 转置矩阵
    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < rows; j++) {
            result[i][j] = matrix[j][i];
        }
    }
}


//卡尔曼增益 Kk = (Pk_ * H')/(H * Pk_ * H' + R)
void Kkcaculate(double pk_[][3], double Kk[3][1])
{
    double H[3] = { 0,0,1 }, H_[3][1] = { {0},{0},{1} };
    double midv1[3][1], midv2[3][1], midv3;
    matrix3v1(pk_, H_, midv1, 3, 3, 1);//
    matrix3v1(pk_, H_, midv2, 3, 3, 1);//
    midv3 = midv2[0][0] * H[0] + midv2[1][0] * H[1] + midv2[2][0] * H[2];
    Kk[0][0] = midv1[0][0] / (midv3 + R);
    Kk[1][0] = midv1[1][0] / (midv3 + R);
    Kk[2][0] = midv1[2][0] / (midv3 + R);

}


//修正估计 x = x_ + Kk * (temp_acc(i)-H*x_);
void xcaclate(double x_[][1], double Kk[][1], double x[][1], double* temp_acc)
{
    double H[3] = { 0,0,1 };
    double mid1, mid2[3];
    mid1 = x_[0][0] * H[0] + x_[1][0] * H[1] + x_[2][0] * H[2];
    mid2[0] = Kk[0][0] * (*temp_acc - mid1);
    mid2[1] = Kk[1][0] * (*temp_acc - mid1);
    mid2[2] = Kk[2][0] * (*temp_acc - mid1);

    x[0][0] = x_[0][0] + mid2[0];
    x[1][0] = x_[1][0] + mid2[1];
    x[2][0] = x_[2][0] + mid2[2];
}

//更新协方差矩阵 Pk = (diag([1 1 1 ]) - Kk * H) * Pk_;
void pkcaculate(double Kk[][1], double pk_[][3], double pk[][3])
{
    double H[3] = { 0,0,1 };
    double mid1[3][3];//存储 Kk * H
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
            {
                mid1[i][j] = 1 - Kk[i][0] * H[j];
            }
            else {
                mid1[i][j] = -Kk[i][0] * H[j];
            }


        }
    }
    matrix3v3(mid1, pk_, pk, 3, 3, 3);

}
/*


 double x[3][1]={{0},{0},{0}};//分别保存位移，速度，加速度
 double pk[3][3]={{1,0,0},{0,1,0},{0,0,1}};//协方差矩阵 Pk

*/
void KalmanFilter(double x[][1], double pk[][3], double* temp_acc)
{
    double dt = 0.004;//采样间隔
    double x_[3][1];

    //状态方程矩阵，换问题就更换矩阵就行
    // double A[3][3]={{1,dt,0.5*dt*dt} , {0,1,dt},{0,0,1}};
    double A[3][3] = { {1,dt,0} , {0,1,dt},{0,0,1} };

    double pk_[3][3], A_[3][3], mid[3][3];
    double Kk[3][1];

    //计算状态方程（先验估计）x_ = A * x ;
    matrix3v1(A, x, x_, 3, 3, 1);

    //先验估计协方差 Pk_ = A * Pk * A' + Q
    transm3x3(A, A_, 3, 3);
    matrix3v3(A, pk, mid, 3, 3, 3);
    matrix3v3(mid, A_, pk_, 3, 3, 3);
    pk_[0][0] += Q;
    pk_[1][1] += Q;
    pk_[2][2] += Q;

    //卡尔曼增益 Kk = (Pk_ * H')/(H * Pk_ * H' + R)
    Kkcaculate(pk_, Kk);

    //修正估计 x = x_ + Kk * (temp_acc(i)-H*x_);
    xcaclate(x_, Kk, x, temp_acc);

    //更新协方差矩阵 Pk = (diag([1 1 1 ]) - Kk * H) * Pk_;
    pkcaculate(Kk, pk_, pk);
}

int main()
{

    FILE* r = fopen("matt.txt", "r");
    //存储卡尔曼滤波后值
    double x[3][1] = { {0},{0},{0} };//分别保存位移，速度，加速度
    double pk[3][3] = { {1,0,0},{0,1,0},{0,0,1} };//协方差矩阵 Pk

    double position[1713],speed[1713],accy[1713];

    //存储加速度原始值
    double acc[1713];

    for (int i = 0; i < 1713; i++)
    {
        fscanf(r, "%lf", &acc[i]);
        
        printf("%lf\n",acc[i]);
    }
    fclose(r);
    


    for (int i = 0; i < 1713; i++)
    {
        
        KalmanFilter( x, pk,&acc[i]);
        position[i]=x[0][0];
        speed[i]=x[1][0];
        accy[i]=x[2][0];
    }
    

    //将结果文件写入文件
    FILE* p = fopen("position.txt", "w");
    FILE* s = fopen("speed.txt", "w");
    FILE* a = fopen("accy.txt", "w");
    for (int i = 0; i < 1713; i++)
    {
        fprintf(p,"%lf\n",position[i]);
        fprintf(s,"%lf\n",speed[i]);
        fprintf(a,"%lf\n",accy[i]);
    }
    fclose(p);fclose(s);fclose(a);
    
   //getchar();
   return 0;
}

