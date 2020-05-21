#include <math.h>
#include <stdio.h>

#define n 7
#define IT_MAX 12

void multiplyMatrices(float A[n][n], float B[n][n], float Res[n][n])
{
	int i, j, k;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			Res[i][j] = 0;
		}
	}

	for(i = 0; i<n; ++i)
	{
		for(j = 0; j<n; ++j)
		{
			for(k=0; k<n; ++k)
			{
				Res[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}


float scalar(float *A, float *B)
{
    float res = 0;
    for (int i = 0; i < n; i++)
        res += A[i] * B[i];
    return res;
}
void matrix_x_vector( float v[n], float M[n][n], float res[n])
{
	for (int i=0; i<n; i++)
    {
        res[i] = 0;
        for (int j=0; j<n; j++)
            res[i] += M[i][j] * v[j];
    }
}
void hotteling(float W[n][n], float lambda, float x_1[n])
{
    float tmp[n][n];
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
        {
            tmp[i][j] = x_1[i] * x_1[j] * lambda;
            W[i][j] = W[i][j] - tmp[i][j];
        }
}
float normalize(float x[n])
{
    return sqrt(scalar(x,x));
}
void transpose(float A[n][n], float A_T[n][n])
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j <n; j++)
            A_T[i][j] = A[j][i];

}

int main()
{
    float A[n][n];
    float W[n][n];
    float D[n][n];
    float result[n][n];
    float x_1[n];
    float x_2[n];

    float lambda[n];
    int i,j;
    int k;


    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
          A[i][j] = 1.0/(sqrt(2 + fabs(i-j)));
          W[i][j] = A[i][j];
        }

    float lambda2;

    for (k=0; k<n; k++)
    {
        for (int i=0; i<n; i++)
            x_1[i] = 1.0;

        for (int i=1; i<=IT_MAX; i++)
        {

            matrix_x_vector(x_1, W, x_2);

            lambda2 = scalar(x_2, x_1) / scalar(x_1, x_1);
            if (k == 0)
                printf("lambda[%d] = %f\n", i, lambda2);

            for (j=0; j<n; j++)
                x_1[j] = x_2[j]/(normalize(x_2));

        }

        for (j=0; j<n; j++)
            result[j][k] = x_1[j];
        lambda[k] = lambda2;

        hotteling(W, lambda[k], x_1);

    }
    puts("");
    puts("Wektory x:");
    for (k=0; k<n; k++)
    {
        printf("x[%d]:\n", k);
        for (j=0; j<n; j++)
            printf("%15.5e", result[j][k]);
        puts("");
    }

    float result2[n][n];
    transpose(result,result2);

    float result3[n][n];
    multiplyMatrices(result2, A, result3);
    multiplyMatrices(result3, result, D);
    puts("Macierz D:");
    for (i=0; i<n; i++)
    {
        printf("  |");
        for (j=0; j<n; j++)
            printf("%15.5e ", D[i][j]);
        printf(" |");
        puts("");
    }

    return 0;
}
