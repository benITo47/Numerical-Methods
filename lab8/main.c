#include <stdio.h>
#include <math.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/gaussj.c"

float f1(float x)
{
    return 1. / (1. + x*x); 
}

float f2(float x)
{
    return cos(2*x);
}

void wyznacz_m(float *x, float *y, float *m, int n);
float wyzSx(float *x, float *y, float *m, int n, float cur_x);
void interpolacja(FILE * fout, float f(float), int n, float x_min, float x_max, float alfa, float beta);
void wyznacz_pochodne(FILE * fout, float f(float), int n, float x_min, float x_max, float alfa, float beta);

int main()
{
    float alfa = 0., beta = 0.;
    float x_min = -5., x_max = 5., step = 0.01;

    FILE * fout_f1 = fopen("f1.dat", "w");
    FILE * fout_f2 = fopen("f2.dat", "w");
    FILE * fout_pochodne = fopen("pochodne.dat", "w");

    interpolacja(fout_f1, f1, 5, x_min, x_max, alfa, beta); // n = 5 
    interpolacja(fout_f1, f1, 8, x_min, x_max, alfa, beta); // n = 8 
    interpolacja(fout_f1, f1, 21, x_min, x_max, alfa, beta); // n = 21

    interpolacja(fout_f2, f2, 5, x_min, x_max, alfa, beta); // n = 5 
    interpolacja(fout_f2, f2, 8, x_min, x_max, alfa, beta); // n = 8 
    interpolacja(fout_f2, f2, 21, x_min, x_max, alfa, beta); // n = 21

    wyznacz_pochodne(fout_pochodne, f1, 10, x_min, x_max, alfa, beta); //pochodne dla f1

    fclose(fout_f1);
    fclose(fout_f2);
    fclose(fout_pochodne);
}

void interpolacja(FILE * fout, float f(float), int n, float x_min, float x_max, float alfa, float beta)
{
    float * x =  vector(1, n);
    float * y =  vector(1, n);
    float * m =  vector(1, n);

    float step = (x_max - x_min) / (n - 1); //krok miedzy wezlami 

    for (int i = 1; i <= n; i++)
        x[i] = y[i] = m[i] = 0.;

    for (int i = 1; i <= n; i++)
    {
        x[i] = x_min + (i-1) * step;
        y[i] = f(x[i]);
    }
    
    wyznacz_m(x, y, m, n);

    for (float cur_x = x_min; cur_x <= x_max; cur_x += 0.01)
        fprintf(fout, "%g\t%g\n", cur_x, wyzSx(x, y, m, n, cur_x));
    
    fprintf(fout, "\n\n");

    free_vector(x, 1, n);
    free_vector(y, 1, n);
    free_vector(m, 1, n);
}


void wyznacz_pochodne(FILE * fout, float f(float), int n, float x_min, float x_max, float alfa, float beta)
{
    float* x = vector(1,n);
    float* y = vector(1,n);
    float* m = vector(1,n);

    float step = (x_max - x_min)/(n-1);

    for(int i = 1; i <= n; i++)
    {
        x[i] = y[i] = m[i] = 0.0; 
    }

    for(int i = 1; i <= n; i++)
    {
        x[i] = x_min + (i-1) * step; 
        y[i] = f(x[i]);
    }

    wyznacz_m(x,y,m,n);// wyznaczamy wartosci drugich pochodnych w wezlach - numerycznie; 
    float dx = 0.01; 

    for(int i = 1; i <=n ; i++ )
        fprintf(fout, "%g\t%g\t%g\n", x[i], m[i], (f(x[i] - dx) - 2*f(x[i]) + f(x[i] + dx)) / (dx*dx)); // x, pochodna^2 numeryczna, pochodna^2 analityczna; 

    free_vector(x, 1, n);
    free_vector(y, 1, n);
    free_vector(m, 1, n);
}
void wyznacz_m(float* x, float*y, float*m, int n)
{
    float** d = matrix(1,n,1,n);
    float** A = matrix(1,n,1,n); 
    float* h = vector(1,n); 
    
    for(int row = 1; row <=n; row++)
    {
        for(int col = 1; col <=n; col++)
        {
            A[row][col] = d[row][col] = 0.0; 
        }
        h[row] = 0.0; 
    }

    for(int i = 2 ; i <=n; i++)
    {
        h[i] = x[i] - x[i-1];
    }

    for(int i = 2; i < n ; i++)
    {
        A[i][i+1] = h[i+1]/(h[i] + h[i+1]); //gorna przekatna 
        A[i][i] = 2.0; //glowna przekatna 
        A[i][i-1] = 1.0 - A[i][i+1]; // dolna przekatna 
    }
    A[1][1] = A[n][n] = 1.0; 
    for(int i = 2; i <n ; i++)
    {
        d[i][1] = 6/ (h[i] + h[i+1]) * ((y[i+1] - y[i])/h[i+1] - (y[i] - y[i-1])/h[i]);
    }

    gaussj(A,n,d,1);

    for(int i = 1 ; i <=n; i++)
    {
        m[i] = d[i][1];
    }
    free_matrix(d, 1, n, 1, n);
    free_vector(h, 1, n);
    free_matrix(A, 1, n, 1, n);

}


float wyzSx(float *x, float *y, float *m, int n, float cur_x)
{
    float* h = vector(1,n); 

    h[1] = 0.0;
    for(int i = 2; i <=n; i++)
    {
        h[i] = x[i] - x[i-1]; //wzor 5 
    }

    int range_number = 0; 

    for(int i = 2; i <= n ; i++ )
    {
        if(cur_x <= x[i])
        {
            range_number = i - 1; 
            break; 
        }
    }

    int i1 = range_number;
    int i = i1 + 1; 

    //wzor 8 
        float s = m[i1] * pow(x[i] - cur_x, 3) / (6. * h[i]) 
             + m[i] * pow(cur_x - x[i1], 3) / (6. * h[i]) 
             + (cur_x - x[i1]) * ((y[i] - y[i1]) / h[i] - h[i] / 6. * (m[i] - m[i1])) 
             + y[i1] - m[i1] * h[i] * h[i] / 6.;

    
    free_vector(h, 1, n);

    return s;
}