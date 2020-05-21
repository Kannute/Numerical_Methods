#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**** Funkcja (2) ****/
double g_x(double x, double y)
{
    return x*x - 4*x + y*y - 4*y + x*y;
}
/**** Rozniczka funkcji (2) po x ****/
double g_deriv_x(double x, double y)
{
    double dx = 0.01;
    double dy = 0.01;
    return ( g_x(x + dx, y) - g_x(x - dx, y) ) / (2.0*dx);
}
/**** Rozniczka funkcji (2) po y ****/
double g_deriv_y( double x, double y)
{
    double dx = 0.01;
    double dy = 0.01;
    return ( g_x(x, y + dy) - g_x(x, y - dy) ) / (2.0*dy);
}
/**** Gradient ****/
void gradient(double *g, double *r_p)
{
    g[0] = g_deriv_x(r_p[0], r_p[1]);
    g[1] = g_deriv_y(r_p[0], r_p[1]);
}

/**** Odwracanie macierzy A ze wzoru****/
void reverse_matrix(double **A)
{
    double a[2][2];
    a[0][0] = A[0][0];
    a[0][1] = A[0][1];
    a[1][0] = A[0][1];
    a[1][1] = A[1][1];

    double tmp = 1.0 / (a[0][0] * a[1][1]  -  a[0][1] * a[1][0]);
    A[0][0] = a[1][1] * tmp;
    A[0][1] = -a[0][1]  * tmp;
    A[1][0] = -a[1][0] * tmp;
    A[1][1] = a[0][0] * tmp;
}
/**** Iloczyn macierz x wektor ****/
void vector_x_matrix(double **A, double *b, double *result)
{
    for (int i=0; i<2; i++)
    {
        result[i] = 0;
        for (int j=0; j<2; j++)
            result[i] += A[i][j] * b[j];
    }
}


/**** Wczytanie nowej pozycji wed³ug wzoru (8) ****/
void new_position(double *r_n, double *r0, double **A, double *g, double omega)
{
    vector_x_matrix(A, g, r_n);
    r_n[0] *= omega;
    r_n[1] *= omega;

    r_n[0] = r0[0] - r_n[0];
    r_n[1] = r0[1] - r_n[1];
}

/**** Implementacja metody Newtona ****/
void newton(double *r0, double omega, double **A)
{
    double r_n[2];
    double g[2];
    int count = 0;
    printf("Dla r0 = [%3.f, %3.f]\n", r0[0], r0[1]);
    //Szukanie dopóty, dopóki metoda nie znajdzie minimum
    while(1)
    {
        count++;
        gradient(g, r0);
        new_position(r_n, r0, A, g, omega);


        //Warunek zakoñczenia
        if (sqrt( pow(r_n[0] - r0[0], 2) + pow(r_n[1] - r0[1], 2) ) < 1e-6)
            break;

        r0[0] = r_n[0];
        r0[1] = r_n[1];
    }
    printf("Ilosc iteracji %d\n", count);
    printf("Polozenie: x:%g, y:%g\n\n", r_n[0], r_n[1]);


}


int main()
{

   /****** 2.1 WYKRES KONTUROWY ******/

    FILE *fptr;
    fptr = fopen("mapa.dat", "a");

    FILE *fptr2;
    fptr2 = fopen("kontur.dat", "a");

    if(fptr == NULL || fptr2 == NULL)
    {
      printf("Error whilst opening a file!\n");
      exit(1);
    }

    //Wpisywanie polozen do pliku
    for (double x = -10.0; x<=10.0; x+= 0.1)
        for (double y = -10.0; y<=10.0; y+=0.1)
        {
            fprintf(fptr, "%g\n", g_x(x,y));
        }


    /**** 2.2 OBLICZANIE MINIMUM ZA POMOCA MACIERZY ODWROTNEJ  ****/
    printf("2.2:\n");

    double **A;
    A = (double **)malloc(2 * sizeof(double *));
    for (int i=0; i<2; i++)
        A[i] = (double *)malloc(2 * sizeof(double));

    A[0][0] = 2.0;
    A[0][1] = 1.0;
    A[1][0] = 1.0;
    A[1][1] = 2.0;

    //Wektor wyrazow wolnych
    double b[2];
    b[0] = -4.0;
    b[1] = -4.0;

    //Wektor wynikowy
    double r[2];

    //Odwrocenie macierzy
    reverse_matrix(A);
    //Obliczenie iloczynu A * b i wpisanie wyniku do wektora r
    vector_x_matrix(A, b, r);

    //We wzorze przed iloczynem stoi -1, wiec mnoze wektor wynikowy
    r[0] = -1.0 * r[0];
    r[1] = -1.0 * r[1];

    printf("Polozenie: x:%g, y:%g\n\n", r[0], r[1]);
	fprintf(fptr2, "%g %g\n\n", r[0], r[1]);
    /**** 2.3 PRZEPIS ITERACYJNY DO ZNALEZIENIA MINIMUM ****/
    printf("2.3:\n\n");

    //Zmienny punkt pocz¹tkowy
    double r0[2];
    //Stala omega
    double omega = 1.0;
    r0[0] = 0;
    r0[1] = 0;

    newton(r0, omega, A);

    r0[0] = 10.0;
    r0[1] = -10.0;
    newton(r0, omega, A);

    r0[0] = 100.0;
    r0[1] = 100.0;
    newton(r0, omega, A);

    r0[0] = 500.0;
    r0[1] = 500.0;
    newton(r0, omega, A);

    /**** 2.4 WPROWADZENIE WAGI < 1 DO PRZEPISU ITERACYJNEGO ****/
    printf("2.4:\n");

    //Staly punkt poczatkowy
    r0[0] = 10.0;
    r0[1] = 10.0;

    //Zmienna omega
    omega = 0.1;
    printf("Dla omega = %.3f\n", omega);
    newton(r0, omega, A);

    r0[0] = 10.0;
    r0[1] = 10.0;

    omega = 0.4;
    printf("Dla omega = %.3f\n", omega);
    newton(r0, omega, A);

    r0[0] = 10.0;
    r0[1] = 10.0;

    omega = 0.7;
    printf("Dla omega = %.3f\n", omega);
    newton(r0, omega, A);


    fclose(fptr);
    fclose(fptr2);
    free(A);
}

