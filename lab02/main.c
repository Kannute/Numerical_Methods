#include <stdio.h>
#include <stdlib.h>

//#include "nrutil.h"
#include "nrutil.h"
#include "nrutil.c"
#include "ludcmp.c"
#include "lubksb.c"
#define N 4

//funkcja liczaca norme macierzy
float max(float **A){
    float maxx = A[1][1];
    for(int i =1; i<= N; i++){
        for(int j=1; j<=N; j++){
           if(maxx<A[i][j])
                maxx = A[i][j];
}

}
return maxx;
}


int main(){
    
    
    //Alokacja macierzy, zmiennych i wektora
    float **A, **A2;
    A = matrix(1, N, 1, N);
    A2 = matrix(1, N, 1, N);
    float d;
    int *indx=ivector(1,N);
    int i,j;
    
    //Wypelniam macierz A i macierz pomocniczaa wartosciami
    for(i=1; i<=N; i++){
        for(j=1; j<=N; j++){
            A[i][j] = ((float)1/(i+j));
            A2[i][j] = ((float)1/(i+j));
    }}
    //uzwywam procedury
    ludcmp(A,N,indx,&d);
    //inicjalizuje wyznacznik i wypisuje elementy diagonalne
    float det=1.0;
    puts("Elementy diagonalne:");
    for(i=1; i<=N; i++){
        j=i;
        printf("%f\n", A[i][j]);
        det *= A[i][j]; 
    }
    det *= d;
    //wypisuje wyznacznik
    printf("wyznacznik A = %g\n", det);
    //inicjalizuje i przypisuje wartosci kolejnym N wektorom
    float *b1=vector(1,N);
     b1[1] =1.0;
     b1[2] =0.0;
     b1[3] =0.0;
     b1[4] =0.0;

    float *b2=vector(1,N);
     b2[2] =1.0;
     b2[1] =0.0;
     b2[3] =0.0;
     b2[4] =0.0;

    float *b3=vector(1,N);
     b3[3] =1.0;
     b3[1] =0.0;
     b3[2] =0.0;
     b3[4] =0.0;

    float *b4=vector(1,N);     
     b4[1] =0.0;
     b4[2] =0.0;
     b4[3] =0.0;
     b4[4] =1.0;
    //inicjalizuje macierz odwrtotna
    float **B = matrix(1,N,1,N);
    //uzywam procedury aby rozwiazac uklad
    lubksb(A,N,indx,b1);
    lubksb(A,N,indx,b2);
    lubksb(A,N,indx,b3);
    lubksb(A,N,indx,b4);
    //wypelniam macierz odwrotna wartosciami
    for(i=1;i<=N;i++){
        B[i][1] = b1[i];
        B[i][2] = b2[i];
        B[i][3] = b3[i];
        B[i][4] = b4[i];
    }
    puts("elementy macierzy odwrotnej:");
    for(i=1; i<=N; i++){       
            printf("%f  %f  %f  %f\n", b1[i],b2[i],b3[i],b4[i] );
    }

    int k;
    float **C=matrix(1,N,1,N);
    for(i=1;i<=N;i++){
        for(j=1;j<=N;j++){
            C[i][j]  =0;
            for(k=1;k<=N;k++)                
                C[i][j] += A2[i][k]*B[k][j];      
        }
    }
    puts("iloczyn: ");

    for(i=1;i<=N;i++){
        for(j=1;j<=N;j++){
            printf("%g\n", C[i][j]);
    }}
    float n1, n2;
    n1 = max(B);
    n2 = max(A2);
    float wsk = n1*n2;
    printf("norma macierzy: : %lf\n", n2 );
    printf("norma macierzy odwrotnej : %lf\n", n1 );
    printf("wskaznik uwarunkowania : %lf\n", wsk );

return 0;
}
