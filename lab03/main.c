#include <stdio.h>
#include <math.h>
#include <time.h>
#define max(X,Y) ((X) > (Y) ? (X) : (Y))
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define abs(X) ((X) > 0 ? (X) : -(X))


//iloczyn skalarny wektorow
double ilWektorow(int length, double *x, double *y) {
	int i = 0;
	double sum = 0.0;
	for(i = 0; i < length; i++)
		sum += x[i] * y[i];

	return sum;
}

//iloczyn Macierz x Wektor
void ilMacierz(int m, int n, double tab[n][n], double *x, double *y) {
	int i = 0, j = 0;

	for(i = 0; i < n; i++) {
		y[i] = 0.0;
		for (j = max(0, i - m); j <= min(i + m, n - 1); j++)
			y[i] += (tab[i][j] * x[j]);

	}
}
//iloczyn Wektor x skalar
void ilSkalar(int n, double *x, double a, double *r) {
	int i = 0;
	for(i = 0; i < n; i++)
		r[i] = (x[i] * a);

}
//odejmowanie wektorow
void odejmowanie(int n, double *x, double *y, double *r) {
	int i = 0;
	for(i = 0; i < n; i++)
		r[i] = x[i] - y[i];

}
//dodawanie wektorow
void dodawanie(int n, double *x, double *y, double *r) {
	int i = 0;
	for(i = 0; i < n; i++)
        r[i] = (x[i] + y[i]);

}

int main() {

	int i = 0, j = 0, k = 0;

	int n = 1000;
	int m = 5;

	double alpha;

	double A[n][n];

	double b[n];
	double r[n];
	double x[n];
	double tmp[n];

    time_t t1,t2;
    double _time;

    //a) dla x == 0
	for (i = 0; i < n; i++)
    {
		x[i] = 0.0;
		b[i] = i+1;

		for(j = 0; j < n; j++)
        {
			if (abs(i - j) <= m)
				A[i][j] = 1.0 / (double)(1 + abs(i - j));
			else
				A[i][j] = 0.0;
		}
	}

	k = 0;
	time(&t1);
	do {

        //tmp =A*x
		ilMacierz(m, n, A, x, tmp);
		//r = b - tmp
		odejmowanie(n, b, tmp, r);

        //tmp = A*r
		ilMacierz(m, n, A, r, tmp);
		//alpha = (r*r)/r*tmp
		alpha = ilWektorow(n, r, r) / ilWektorow(n, r, tmp);

        //tmp = alpha*r
		ilSkalar(n, r, alpha, tmp);
		//x = x + tmp
		dodawanie(n, x, tmp, x);

//wypisanie kroku, wartosc normy euklidesowej wektora reszt, alpha oraz wartosc normy euklidesowej wektora rozwiazan
		printf("%d ", k);
		printf("%f ", sqrt(ilWektorow(n, r, r)));
        printf("%f ", alpha);
        printf("%f\n", sqrt(ilWektorow(n, x, x)));

		k++;
	} while( sqrt( ilWektorow(n, r, r) ) > 1e-6);

	//wypisanie czasu
    time(&t2);
    _time=difftime(t1,t2);
    printf("czas potrzebny na wykonanie : %f\n", _time);

	//b) dla x == 1
	for (i = 0; i < n; i++)
		x[i] = 1.0;

	//wszystko analogicznie do podpunktu a)
    k=0;
	do {

		ilMacierz(m, n, A, x, tmp);
		odejmowanie(n, b, tmp, r);

		ilMacierz(m, n, A, r, tmp);
		alpha = ilWektorow(n, r, r) / ilWektorow(n, r, tmp);

		ilSkalar(n, r, alpha, tmp);
		dodawanie(n, x, tmp, x);

        printf("%d ", k);
		printf("%f ", sqrt(ilWektorow(n, r, r)));
        printf("%f ", alpha);
        printf("%f\n", sqrt(ilWektorow(n, x, x)));

		k++;
	} while( sqrt( ilWektorow(n, r, r) ) > 1e-6);



	return 0;
}
