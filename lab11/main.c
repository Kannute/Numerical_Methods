#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

/** Sygnal okresowy niezaszumiony **/
double y_0(int i, int N)
{
    double omega = 2*(2*M_PI)/N;
    return sin(omega*i) + sin(2*omega*i) + sin(3*omega*i); 
}

/** Zmienna losowa imitujaca szum **/
double generate_delta()
{
    return 2*(rand()/(RAND_MAX + 1.0)) - 1;
}

/** Sygnal zaszumiony **/
double y_i(int i, int N)
{
    return y_0(i, N) + generate_delta();
}

/** Modul z liczby zespolonej **/
double complex_module(double a, double b)
{
    return sqrt(pow(a,2) + pow(b,2));
}

/** Wykonanie wszystkiego dla zadanego k **/
void step(int k)
{
    
    int N = pow(2,k);


    //Alokacja wektora szumu y
    gsl_complex_packed_array y = malloc(2*N * sizeof(double));

    //Wypelnienie go odpowiednimi wartosciami
    for(int i=0;i<N;i++)
    {
        y[2*i] = y_i(i, N);
        y[2*i +1] = 0.0;
    }

    //Wypisanie
    for(int i=0;i<N;i++)
        printf("%d %lf\n", i, y[2*i]);
    printf("\n\n");

    if(k==8)
    {
        FILE *fp8;
        fp8 = fopen("y8.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp8, "%d %lf\n", i, y[2*i]);
        fprintf(fp8, "\n\n");
        fclose(fp8);
    }
    if(k==10)
    {
        FILE *fp10;
        fp10 = fopen("y10.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp10, "%d %lf\n", i, y[2*i]);
        fprintf(fp10, "\n\n");
        fclose(fp10);
    }
    if(k==12)
    {
        FILE *fp12;
        fp12 = fopen("y12.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp12, "%d %lf\n", i, y[2*i]);
        fprintf(fp12, "\n\n");
        fclose(fp12);
    }

    
    
    //FFT
    gsl_fft_complex_radix2_forward(y, 1, N);

    if(k==8)
    {
        FILE *fp;
        fp = fopen("fft8.dat", "a");
        for (int i = 0; i < N; i++)
            fprintf(fp, "%d %g %g\n", i, y[2 * i], y[2 * i + 1]);
        fclose(fp);
    }
    double max = complex_module(y[0], y[1]);
    double tmp;
    
    //Szukanie najwiekszej wartosci
    for(int i=2; i<2*N; i+=2)
    {
        tmp = complex_module(y[i], y[i+1]);
        if(max < tmp)
            max = tmp;
    }

    printf("max = %lf\n\n\n", max);
    max*=0.5;

    //Dyskryminacja
    for(int i=0; i<2*N; i+=2)
    {
        tmp = complex_module(y[i], y[i+1]);
        if(tmp < max)
        {
            y[i] = 0.0;
            y[i+1] = 0.0;
        }

    }
    //FFT part 2
    gsl_fft_complex_radix2_backward(y, 1, N);

    //Podzielenie wynikow przez N
    for(int i=0;i<N;i++)
    {
        y[2*i] = y[2*i] / ( (double)N );
        printf("%d %lf\n", i, y[2*i]);
    }

    if(k==8)
    {
        FILE *fp8;
        fp8 = fopen("y8.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp8, "%d %lf\n", i, y[2*i]);
        fprintf(fp8, "\n\n");
        fclose(fp8);
    }
    if(k==10)
    {
        FILE *fp10;
        fp10 = fopen("y10.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp10, "%d %lf\n", i, y[2*i]);
        fprintf(fp10, "\n\n");
        fclose(fp10);
    }
    if(k==12)
    {
        FILE *fp12;
        fp12 = fopen("y12.dat", "a");
        for(int i=0;i<N;i++)
            fprintf(fp12, "%d %lf\n", i, y[2*i]);
        fprintf(fp12, "\n\n");
        fclose(fp12);
    }
    printf("\n\n");
    
    //Zwalnianie pamieci
    free(y);
   
    
}

int main()
{
    //Dla k==8
    step(8);

    //Dla k==10
    step(10);

    //Dla k==12
    step(12);

    return 0;
}

