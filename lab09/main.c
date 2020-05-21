#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

double x0 = 2.0;
double sigma = 4.0;
double xmin = -10.0;
double xmax = 14.0;

//Wyznaczanie wektora polozen rownoodleglych wezlow
void x_w(double N, double *xw)
{
    double step = (xmax - xmin) / (N - 1);
    for (int i=0; i<N; i++)
        xw[i] = xmin + i * step;
    
}

//Wyznaczanie pojedynczej wartosci wektora wartosci funkcji g(x)
double g_x_i(double x_w_i, double alfa)
{
    double a0, a1, a2;
    a0 = -pow(x0, 2) / (2*pow(sigma, 2));
    a1 = x0 / pow(sigma, 2);
    a2 = -1 / (2*pow(sigma,2));

    double U = rand()/(RAND_MAX + 1.0);

    double delta = alfa * (U - 0.5);

    return exp(a0 + a1*x_w_i + a2*pow(x_w_i,2)) * (1 + delta);
}

//Wyznaczanie wektora wartosci funkcji g(x) na wezlach g(xj)
void g_w(double N, double *gw, double *xw, double alfa)
{
    for (int i=0; i<N; i++)
        gw[i] = g_x_i(xw[i], alfa);
}

//Wyznaczenie wketora wartosci funkcji f(x) na podstawie poprzedniego wektora g(w)
double f_w(double N, double *gw, double *fw)
{
    for (int i=0; i<N; i++)
        fw[i] = log(gw[i]);
}

//wektor wyrazow wolnych r o m elementach, korzystajac ze wzoru (8)
void make_r(gsl_vector *r, double *xw, double *fw, double m, double N)
{
    double step;
    for (int k=0; k<m; k++)
    {
        step = 0;

        for (int j=0; j<N; j++)
            step += fw[j] * pow(xw[j], k);

        gsl_vector_set(r, k, step);
    }
}

//Wyznaczanie macierzy G o rozmiarze mxm wypelniana zgodnie ze wzorem (9)
void make_G(gsl_matrix *G, double *xw, double m, double N)
{
    double step;
    for (int i=0; i<m; i++)
    {
        for (int k=0; k<m; k++)
        {
            step = 0;

            for (int j=0; j<N; j++)
                step += pow(xw[j], i+k);

            gsl_matrix_set(G, i, k, step);
        }
    }
}

//Wartosc funkcji aproksymujacej
double G_x_i(double x, double *b)
{
    return exp(b[0] + b[1] * x + b[2] * pow(x,2) + b[3] * pow(x,3));
}

//Wyznaczanie wektora b
void make_b(gsl_vector *r, gsl_matrix *G, double *b, double m)
{
   //Householder
    gsl_linalg_HH_svx(G, r);

    for (int i=0; i<m; i++)
        b[i] = gsl_vector_get(r, i);
}


//Wypisywanie do pliku pkt.dat
void pkt_dat(double *xw, double *gw,int N)
{
    FILE *f;
    f = fopen("pkt.dat", "a");
    if(!f)
    {
        puts("Blad");
        exit(-1);
    }
    for (int i=0; i<N; i++)
        fprintf(f, "%-15g %-15g\n", xw[i], gw[i]);
    fprintf(f, "\n\n");
    fclose(f);
}

//Wypisywanie do pliku G.dat
void G_dat(double step, double *b)
{
    FILE *f;
    f = fopen("G.dat", "a");
    if(!f)
    {
        puts("Blad");
        exit(-1);
    }
    for(double x=xmin; x<xmax; x += step)
        fprintf(f, "%-15g %-15g\n", x, G_x_i(x, b));
    fprintf(f, "\n\n");
    fclose(f);    
}

void full(double N, double alfa)
{
    //Alokacja wektora wezlow  i wypelnienie go wartosciami
    double *xw = (double*)malloc(N*sizeof(double));
    x_w(N, xw);
    
    //Alokacja wektora g(w) i wypelnienie go wartosciami
    double *gw = (double*)malloc(N*sizeof(double));
    g_w(N, gw, xw, alfa);

    //Alokacja wektora fw i wypelnienie go wartosciami
    double *fw = (double*)malloc(N*sizeof(double));
    f_w(N, gw, fw);
    
    //Wpisanie wartosci punktowych do pliku
    pkt_dat(xw, gw, N);

    double m = 4;
    //Wektor reszt
    gsl_vector *r = gsl_vector_calloc(m);
    make_r(r, xw, fw, m, N);
    
    //Macierz G
    gsl_matrix *G = gsl_matrix_calloc(m, m);
    make_G(G, xw, m ,N);

    //Wektor wynikowy
    double *b = (double*)malloc(m*sizeof(double));
    make_b(r, G, b, m);

    //Wpisanie wynikow do pliku
    G_dat(0.1, b);

    free(xw);
    free(gw);
    free(fw);
    free(b);
    gsl_matrix_free(G); 
    gsl_vector_free(r);
}


int main()
{
    //Dla 11 wezlow bez zaburzen
    int N = 11;

    double alfa = 0;
    full(N, alfa);
    
    //przypadek z losowym zaburzeniem dla 11 wezlow
    alfa = 0.5;
    full(N, alfa);
    
    //Dla 101 wezlow
    N = 101;
    full(N, alfa);


    return 0;
}