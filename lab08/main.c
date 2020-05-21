#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//Krance 
#define xmax 5
#define xmin -5

//Prototypy
void wyznacz_M(double *xw, double *yx, gsl_vector *m, int n, double alfa, double beta);
double wyznacz_Sx(double *xw, double *yx, gsl_vector *m, int n, double x);

//Dane funkcje f1(x) oraz f2(x)
double f_1(double x)
{
    return 1.0/(1.0 + pow(x, 2));
}

double f_2(double x)
{
    return cos(2.0*x);
}
double deriv(double x, double dx)
{
    return (f_1(x - dx) - 2*f_1(x) + f_1(x + dx)) / pow(dx, 2);
}
void printVector(gsl_vector *v, int n);
void printVectorAcc(gsl_vector *v, int n);
void printMatrix(gsl_matrix *M, int n);

int main()
{
    /*
    FILE *f1;
    f1 = fopen("f1.dat","w");
    if(f1 == NULL)
        exit(1);
    */
    FILE *fdrv;
/*
    FILE *f2;    
    f2 = fopen("f2.dat","w");
    if(f2 == NULL)
       exit(1);
    */
    fdrv = fopen("pochodne.dat","w");
    if(fdrv == NULL)
    {
        printf("Error!");   
        exit(1);             
    }

    
    //Ustalone wartosci alfa i beta
    double alfa = 0;
    double beta = 0;
    //Liczba wezlow. Poczatkowo rowna 5
    int n = 10;
    //Tablica polozen wezlow
    double xw[n];
    //Tablica wartosci funkcji
    double yx[n];
    //Stala delta
    double deltax = ((double)xmax - (double)xmin) / ((double)n - 1.0);
    int i,f;
    //delta x potrzebna do liczenia pochodnych podwojnych
    double dx = 0.01;

    //Wypelnianie tablcy wedlug wzoru (11)
    for(i=0;i<n;i++)
        xw[i] = (double)xmin + deltax*i;
        
    //Wypelnienie wektora wartosci funkcji 
    for(i=0;i<n;i++)
    {
        //yx[i] = f_2(xw[i]);
        //f = 2;
        yx[i] = f_1(xw[i]);
        f = 1;
    }

    /*printf("\nDla yx[i] = f%d(x[i]).\n" , f);
    for(i=0;i<n;i++)
        printf("| %d | %f | %1.8f |\n" , i, xw[i], yx[i] );
    */
    //Alokacja wektora wypelnionego zerami
    gsl_vector *m = gsl_vector_calloc(n);
    
    wyznacz_M(xw, yx, m, n, alfa, beta);
/*
    for(double x = (double)xmin; x<= (double)xmax; x+= 0.01) 
    { 
        double temp = wyznacz_Sx(xw, yx, m, n, x);
        fprintf(f2,"%f %f\n",x, temp);
        //printf("%1.3f %2.7e\n", x, temp);
    }
    fprintf(f2, "\n\n");*/
    
    if(f == 1 && n == 10)
    {   
        double temp;
        for(i=0;i<n;i++)
        {
            temp = deriv(xw[i] , dx);
            fprintf(fdrv,"%f %f %f\n", xw[i], gsl_vector_get(m, i), temp);
            //printf("%f %f %f\n", xw[i], gsl_vector_get(m, i), temp);
        }
    } 

    gsl_vector_free(m);
    
    //fclose(f2);
    //fclose(f2);
    fclose(fdrv);
    
    return 0;
}

void wyznacz_M(double *xw, double *yx, gsl_vector *m, int n, double alfa, double beta)
{
    //Alokacja macierzy ukladu, wektora wyrazow wolnych oraz zmiennych
    gsl_matrix * A = gsl_matrix_alloc(n, n);
    gsl_vector * d = gsl_vector_alloc(n);

    //Odleglosci miedzywezlowe sa w naszym przypadku stale, wiec h jest stale
    double h_i = ((double)xmax - (double)xmin) / ((double)n - 1.0);
    //co za tym idzie stale i niezalezne od n sa zmienne lambda i mu
    double lambda = h_i / (h_i + h_i);
    double mu = 1.0 - lambda;
    
    //Ustawiam wartosci alfa oraz beta kolejno na poczatku i koncu wektora
    gsl_vector_set(d, 0, alfa);
    gsl_vector_set(d, n-1, beta);

    //Pozostale wartosci ustawiam wedlug wzoru (7)
    double value;
    for(int i = 1; i < n-1; i++)
    {
        value = (6.0 / (h_i + h_i)) * ( (yx[i+1] - yx[i]) / h_i - (yx[i] - yx[i-1]) / h_i);
        gsl_vector_set(d, i, value);
    }
    
    //printf("\nWektor d dla n = %d\n" , n);
    //printVector(d,n);
    
    
    //Ustawiam wszystkie elementy macierzy A na zero
    gsl_matrix_set_zero(A);

    //Ustawiam diagonale wedlug wzoru (7)
    for(int i = 0; i < n; i++)
        gsl_matrix_set(A, i, i, 2);  
        
    //Ustawianie diagonali cont.    
    gsl_matrix_set(A, 0, 0, 1);
    gsl_matrix_set(A, n-1, n-1, 1);

    //Wypelnianie odpowiednimi wartosciami ze wzoru (7) elementow pod i nad diagonala
    for(int i = 1; i < n-1; i++)
    {
        for(int j = 1; j < n-1 ; j++)
        {
            if(i == j)
            {
                gsl_matrix_set(A, i, j-1, mu);
                gsl_matrix_set(A, i, j+1, lambda);
            }
        }   
    }
    
    //printf("\nMacierz A dla n = %d" , n);
    //printMatrix(A,n);

    //transformacja Householder'a
    gsl_linalg_HH_svx( A, d);

    //Przypisanie wektorowi m wartosci wektora wynikowego d
    for(int i=0; i<n;i++)
        gsl_vector_set(m, i, gsl_vector_get(d, i));
    
    //printf("\nWektor m\n");
    //printVectorAcc(m, n);
    //zwalnianie pamieci
    gsl_matrix_free(A); 
    gsl_vector_free(d);

}

double wyznacz_Sx(double *xw, double *yx, gsl_vector *m, int n, double x)
{
    int przedzial;
    //zmienna wynikowa Sx
    double Sx;
    //stala
    double h_i = ((double)xmax - (double)xmin) / ((double)n - 1.0);

    //Szukanie danego przedzialu = i-1
    for(int i=1; i<n; i++)
    {
        if(xw[i-1] <= x && x <= xw[i])
        {
            przedzial = i-1;
            break;
        }
    }

    //Wyznaczanie stalych A_i oraz B_i
    double A_i;
    double B_i;
    
    double tmp1 = gsl_vector_get(m , przedzial + 1);
    double tmp2 = gsl_vector_get(m , przedzial);
    A_i =( (yx[przedzial+1] - yx[przedzial]) / h_i ) - (h_i/6.0) * ( tmp1 - tmp2);    

    tmp1 = gsl_vector_get(m , przedzial);
    B_i = yx[przedzial] - tmp1 * ( (pow(h_i,2)) / 6.0 );

    //Wyznaczanie wzoru (8), stopniowo, zeby sie nie pogubic
    Sx = gsl_vector_get(m , przedzial);
    Sx *= (pow((xw[przedzial+1] - x), 3) / (6.0*h_i));
    Sx += gsl_vector_get(m , przedzial + 1) * (pow((x - xw[przedzial]), 3) / (6.0*h_i));
    //Sx *= pow((x - xw[przedzial]), 3) / 6.0*h_i;
    Sx += A_i * (x - xw[przedzial]);
    //Sx *= (x - xw[przedzial]);
    Sx += B_i;

    return Sx;


}

void printVector(gsl_vector *v, int n)
{
    puts("");
    double temp;
        for(int i=0; i<n ; i++)
        {
            temp =  gsl_vector_get(v , i);
            printf("| %.7f |\n", temp);     
        }

}

void printVectorAcc(gsl_vector *v, int n)
{
    puts("");
    double temp;
        for(int i=0; i<n ; i++)
        {
            temp =  gsl_vector_get(v , i);
            printf("| %.7e |\n", temp);     
        }
}

void printMatrix(gsl_matrix *M, int n)
{
    double temp;
    for(int i=0; i<n ; i++)
        {
            puts("");
            for(int j=0;j<n;j++)
            {
                temp =  gsl_matrix_get(M , i, j);
                //printf(" %2.2f ", temp);
                printf("Value for [%d][%d] = %2.2f\n",i,j, temp); 
            }
        }
    puts("");
}




