#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/gauleg.c"
#include "/home/NR/numerical_recipes.c/gauher.c"
#include "/home/NR/numerical_recipes.c/gaulag.c"
#include "/home/NR/numerical_recipes.c/gammln.c"
#include <stdio.h>
#include <stdlib.h> 

float func1(float x)
{
    return 1.0/(x * sqrt(x*x - 1));
}

float func2(float x)
{
    return log(x)*exp(-x*x);
}

float func3(float x)
{
    return sin(2*x) * exp(-2 * x ); // originaly was sin(2*x)exp(-3 * x )
}

float func3b(float x)
{
    return 1.0/3.0 * sin(2.0/3.0 * x);
}


void calculateFirstTask(FILE* ptr, float a, float b)
{      
    float* x = NULL; 
    float* w = NULL; 
    float c = 0.0;  

    for(int i = 2; i <= 100; i++)
    {
        int n = i; 
        x = vector(1,n);
        w = vector(1,n);
        gauleg(a,b,x,w,n);

        
        c = 0.0;
        for(int i = 1; i <= n; i++)
        {
          c +=  w[i] * func1(x[i]);
        }

        fprintf(ptr, "%d\t%g\n",n, fabs(c - (M_PI/3.0))); 
        free_vector(x,1,n);
        free_vector(w,1,n);
    }
    
    fprintf(ptr,"\n\n\n\n");
}

void calculateSecondTask(FILE* ptr)
{
    
    float* x = NULL; 
    float* w = NULL; 
    float c = 0.0;  
    //Hermite
    for(int i = 2; i <= 100; i+=2)
    {
        int n = i; 
        x = vector(1,n);
        w = vector(1,n);
        gauher(x,w,n);

        
        c = 0.0;
        for(int i = 1; i <= n; i++)
        {
          c +=  w[i] * 0.5 * log(fabs(x[i]));
        }

        fprintf(ptr, "%d\t%g\n",n, fabs(c + 0.8700577)); 
        free_vector(x,1,n);
        free_vector(w,1,n);
    }
    

    fprintf(ptr,"\n\n\n\n");
    //Legendre
    for(int i = 2; i <= 100; i++)
    {
        int n = i; 
        x = vector(1,n);
        w = vector(1,n);
        gauleg(0,5.0,x,w,n);

        
        c = 0.0;
        for(int i = 1; i <= n; i++)
        {
          c +=  w[i] * func2(x[i]);
        }

        fprintf(ptr, "%d\t%g\n",n, fabs(c + 0.8700577)); 
        free_vector(x,1,n);
        free_vector(w,1,n);
    }
    
    fprintf(ptr,"\n\n\n\n");

}

void calculateThirdTask(FILE* ptr)
{
    float* x = NULL; 
    float* w = NULL; 
    float c = 0.0;  
    
    for(int i = 2; i <= 20; i++)
    {
        int n = i; 
        x = vector(1,n);
        w = vector(1,n);

        gaulag(x,w,n,0.0);

        c = 0.0;
        for(int i = 1; i <= n; i++)
        {
          c +=  w[i] * func3(x[i]);
        }

        fprintf(ptr, "%d\t%g\n",n, fabs(c - (2.0/13.0))); 
        free_vector(x,1,n);
        free_vector(w,1,n);
    }
    fprintf(ptr,"\n\n\n\n");

        for(int i = 2; i <= 20; i++)
    {
        int n = i; 
        x = vector(1,n);
        w = vector(1,n);

        gaulag(x,w,n,0.0);

        c = 0.0;
        for(int i = 1; i <= n; i++)
        {
          c +=  w[i] * func3b(x[i]);
        }

        fprintf(ptr, "%d\t%g\n",n, fabs(c - (2.0/13.0))); 
        free_vector(x,1,n);
        free_vector(w,1,n);
    }
    fprintf(ptr,"\n\n\n\n");

}

int main()
{
    FILE* fptr1 = fopen("out.dat", "w");
    calculateFirstTask(fptr1, 1.0,2.0);

    calculateSecondTask(fptr1);

    calculateThirdTask(fptr1);

    fclose(fptr1);
}   