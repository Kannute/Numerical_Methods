#include <stdio.h>
#include <math.h>


//Wielomian
double f_x(double x)
{
    return (x - 1.2) * (x - 2.3) * pow((x - 3.3), 2);
}

//Pochodna wielomianu
double drv_f_x(double x, double delta_x)
{
    return (f_x(x + delta_x) - f_x(x - delta_x)) / (2*delta_x);
}

//Element xk+1 dla metody niemofyfikowanej
double xk_plus_1(double xk, double xk_minus_1)
{
    return xk - ((xk - xk_minus_1) * f_x(xk)) / (f_x(xk) - f_x(xk_minus_1));
}

//Zastepujaca f(x) funkcja u(x)
double u_x(double x, double delta_x)
{
    return f_x(x) / drv_f_x(x , delta_x);
}
//Element xk+1 dla metody modyfikowanej
double xk_plus_1_ver_mod(double xk, double xk_minus_1, double delta_x)
{
    return xk - ((xk - xk_minus_1) * u_x(xk, delta_x )) / (u_x(xk , delta_x) - u_x(xk_minus_1, delta_x));
}


int main()
{
    double x0, x1, x2;
    int k;
    printf("Niemodyfikowana metoda siecznych:\n\n");

    /*Pierwszy przypadek*/
    x0 = 0.9; x1 = 1.0;
    printf("Dla pary x0= %.2f, x1 = %.2f\n", x0, x1);

    for (k = 1; k <= 10000; ++k)
    {
        x2 = xk_plus_1(x1,x0);
        printf("k = %d, x2 = %f, e = %6.7e, f(x2) = %6.7e\n", k, x2, fabs(x2 - x1), f_x(x2));
        if(fabs(x2-x1) < 1e-6)
            break;
        x0 = x1;
        x1 = x2;
    }

    printf("\n\n\n");

    /*Drugi przypadek*/
    x0 = 1.7; x1 = 1.75;
    printf("Dla pary x0= %f, x1 = %f\n", x0, x1);

    for (k = 1; k <= 10000; ++k)
    {
        x2 = xk_plus_1(x1,x0);
        printf("k = %d, x2 = %f, e = %6.7e, f(x2) = %6.7e\n", k, x2, fabs(x2 - x1), f_x(x2));
        if(fabs(x2-x1) < 1e-6)
            break;
        x0 = x1;
        x1 = x2;
    }

    printf("\n\n\n");

    /*Trzeci przypadek*/
    x0 = 3.7; x1 = 3.65;
    printf("Dla pary x0= %f, x1 = %f\n", x0, x1);


    for (k = 1; k <= 10000; ++k)
    {
        x2 = xk_plus_1(x1,x0);
        printf("k = %d, x2 = %f, e = %6.7e, f(x2) = %6.7e\n", k, x2, fabs(x2 - x1), f_x(x2));
        if(fabs(x2-x1) < 1e-6)
            break;
        x0 = x1;
        x1 = x2;
    }

    printf("\n\n\n");



    printf("Modyfikowana metoda siecznych ");

    double delta_x;

    /*Przypadek dla delty rownej 0.1*/
    x0 = 1.7; x1 = 1.75;
    //x0 = 0.9; x1 = 1.0;
    //x0 = 1.7; x1 = 1.75;

    printf("Dla x0 = %.2f, x1 = %.2f\n\n", x0, x1);

    delta_x = 0.1;
    printf("Delta x = %.1f \n", delta_x);

    for (k = 1; k <= 10000; ++k)
    {
        x2 = xk_plus_1_ver_mod(x1,x0, delta_x);
        printf("k = %d, x2 = %f, e = %6.7e, f(x2) = %6.7e\n", k, x2, fabs(x2 - x1), f_x(x2));
        if(fabs(x2-x1) < 1e-6)
            break;
        x0 = x1;
        x1 = x2;
    }

    printf("\n\n\n");

    /*Przypadek dla delty rownej 0.001*/
    x0 = 0.9; x1 = 1.0;
    //x0 = 0.9; x1 = 1.0;
    //x0 = 1.7; x1 = 1.75;

    delta_x = 0.001;
    printf("Delta x = %.3f \n", delta_x);

    for (k = 1; k <= 10000; ++k)
    {
        x2 = xk_plus_1_ver_mod(x1,x0, delta_x);
        printf("k = %d, x2 = %f, e = %6.7e, f(x2) = %6.7e\n", k, x2, fabs(x2 - x1), f_x(x2));
        if(fabs(x2-x1) < 1e-6)
            break;
        x0 = x1;
        x1 = x2;
    }

}
