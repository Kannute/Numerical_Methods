#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X_MAX -5.0
#define X_MIN 5.0
#define PI 3.14

double F_x(double x)
{
    return exp(-pow(x, 2));
}

double Lagrange_element(double x, int n, double *x_m ,double *y_m)
{
    double w_x = 0.0;
    double a;
    int j,k;

    for (j=0; j<=n; j++)
    {
        a = 1.0;
        for (k=0; k<=n; k++)
        {
            if (k!=j)
                a *= ((x - x_m[k]) / (x_m[j]-x_m[k]));
        }
        w_x += y_m[j] * a;
    }
    return w_x;
}

void czebyszew(double *x_m , int n)
{
    for (int m=0; m<=n; m++)
        x_m[m] = (0.5) * ( (X_MAX - X_MIN) * cos(PI * (double)(2*m + 1) / (double)(2*n + 2) ) + (X_MIN + X_MAX) );
}

int main()
{
    int n =20;
    double x_m[n+1];
    double y_m[n+1];
    double W[1001];
    double h;
    int i;


    h = (X_MIN - X_MAX)/n;

    x_m[0] = X_MAX;
    x_m[n] = X_MIN;

    for(i=1; i<n; i++)
        x_m[i] = X_MAX + i*h;

    for (i=0; i<=n; i++)
        y_m[i] = F_x(x_m[i]);

    printf("Rownolegle wiezy interpolacji dla n=%d\n", n);
    puts("| i | x_m | y_m |");

    for (i=0; i<=n; i++)
        printf("| %d | %1.1f | %2.5e |\n", i, x_m[i], y_m[i]);



    czebyszew(x_m,n);

    for (i=0; i<=n; i++)
        y_m[i] = F_x(x_m[i]);

    printf("\n\n");

    printf("Wezly interpolacji - zera wielomianow czebyszewa dla n= %d\n", n);
    for (i=0; i<=n; i++)
        printf("| %d | %1.1f | %2.5e |\n", i, x_m[n-i], y_m[n-i]);


    printf("\n\n");
    puts("Wielomian interpolacyjny Lagrange'a:");

    for (i=0; i<1001; i++)
    {
        W[i] = Lagrange_element(i*0.01 - 5, n, x_m, y_m);
        //printf("W%d(%2.2f) = %2.7e\n", i, i*0.01 -5, W[i]);
    }




    return 0;
}
