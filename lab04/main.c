#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "nr.h"

#define n 5
void MatrixTransposition(float **A)
{
    float **temp;
    temp = matrix(1,n,1,n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <=n; j++)
            temp[i][j] = A[j][i];

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            A[i][j] = temp[i][j];

    free(temp);
}

float** multiplyMatrixes(float **A, float **B)
{
    float **X;
    X = matrix(1, n, 1, n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            X[i][j] = 0.0;

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            for (int k = 1; k <= n; k++)
                X[i][j] += A[i][k] * B[k][j];

    return X;
}

float Scalar(float *A, float *B)
{
    float res = 0;
    for (int i = 1; i <= n; i++)
        res += A[i] * B[i];
    return res;
}



int main(){
    //inicjalizacja wektorow i wypelnienie ich odpowiednimi wartosciami
	float** A = matrix(1,n,1,n);
	float *d = vector(1, n);
	float *e = vector(1,n);

	int i, j;
	for(i=1; i<=n; i++)
		for(j=1; j<=n; j++)
			A[i][j] = sqrt(i+j);

    for(i=1; i<=n; i++)
    {
        d[i] = 0;
        e[i] = 0;
    }



    //redukcja householder'a i wypisanie macierzy przeksztalcenia P
	tred2(A, n, d,e);
	printf("\nMacierz przeksztalcenia P:\n");
	for(i=1; i<=n; i++)
    {
		for(j=1; j<=n; j++)
			printf("%17.7f ", A[i][j]);
		puts("");
	}




    //Inicjalizacja i wypelnienie macierzy jednostkowej Y
	float** Y = matrix(1,n,1,n);

	for(i=1; i<=n; i++)
        for(j=1; j<=n; j++)
            Y[i][j] = 0;
	for(i=1; i<=n; i++)
		Y[i][i] = 1;

    //rozwi¹zanie rownania wlasnego dla rzeczywistej symetrycznej macierzy trójdiagonalnej

	tqli(d, e, n, Y);
   // d[1] = -0.00000024575;
	printf("\nMacierz Y po wykonaniu funkcji tqli. Kolumny przechowuja wektory wlasne T:\n");
	for(i=1; i<=n; i++)
        {
		for(j=1; j<=n; j++)
			printf("%17.7e ", Y[i][j]);
		puts("");
        }

    //wypisanie wartosci wlasnych dla T i A
    printf("\nWartosci wlasne macierzy A\n");
    for(i=1; i<=n; i++)
        printf("%17.7e ", d[i]);
    puts("");


	//inicjalizacja macierzy wektorow wlasnych X. Gdzie X = A*Y
	float** X =  multiplyMatrixes(A,Y);

    printf("\nWektory wlasne macierzy A:\n");
    for(i=1; i<=n; i++){
		for(j=1; j<=n; j++)
			printf("%17.7e ", X[i][j]);
		puts("");
	}

    //inicjalizacja macierzy pomocniczej A_tmp = A
    float** A_tmp = matrix(1,n,1,n);
	for(i=1; i<=n; i++)
		for(j=1; j<=n; j++)
			A_tmp[i][j] = sqrt(i+j);


    float **Z = multiplyMatrixes(A_tmp, X);
    //transpozycja macierzy
    MatrixTransposition(Z);
    MatrixTransposition(X);

    //inicjalizacja wektora wynikowego ilorazu sum
    float *beta =  vector(1, n);

    for(i = 1; i <= n; i++)
        beta[i] = Scalar(Z[i], X[i]) / Scalar(X[i], X[i]);


    printf("\nwartosci beta:\n");
    for(i = 1; i <= n; i++)
        printf("%17.7e", beta[i]);
    puts("");


    free_matrix(A,1,n,1,n);
    free_matrix(Z,1,n,1,n);
    free_matrix(X,1,n,1,n);
    free_matrix(Y,1,n,1,n);
    free_matrix(A_tmp,1,n,1,n);
    free_vector(beta, 1 , n);
    free_vector(d, 1, n);
    free_vector(e,1,n);


	return 0;
}
