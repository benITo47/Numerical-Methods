#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define n 8

double f(double x);
double Simpson(FILE * ptr);
double Milne(FILE * ptr); 
void PrintToScreen(int size, double TAB[][size]);

int main(){
    FILE * ptr = fopen("result.dat","w");

    double resultSimpson = Simpson(ptr);
    fprintf(ptr,"\n\n");
    double resultMilne = Milne(ptr);

    printf("\nMilne: %.10g\n", resultSimpson);
    printf("\nSimpson: %.10g\n", resultSimpson);

    fclose(ptr);
}

double f(double x) {
    return log(x*x*x + 3*x*x + x + 0.1)*sin(18*x);
}

double Simpson(FILE * ptr) {
    double D[n+1][n+1], h;
    int N;     

    fprintf(ptr, "Simpson:\n");

    for(int i=0; i<n+1; i++)
    {
        D[i][0] = 0.;
        N = pow(2,i+1);
        h = 1./(double)N;
        for(int j=0; j<=N/2 - 1; j++)
            D[i][0] += h/3.*(0+f(2*j*h) + 4*f(0+(2*j+1)*h) + f(0+(2*j+2)*h));
        for(int k=1; k<i+1; k++)
            D[i][k] = (pow(4,k)*D[i][k-1] - D[i-1][k-1])/(pow(4,k) - 1);
        fprintf(ptr,"%d\t%15.10g\t%15.10g\n", i, D[i][0], D[i][i]);
    }

    printf("\n\n------------------------Simpson------------------------\n");
    printf("D[w][0] \t D[w][w] \n");
    for(int i =0; i <n ;i++)
    {
            printf("%.10g \t %.10g \n", D[i][0], D[i][i]);
    }
    
    printf("\n\n");

    for(int i = 0; i< n; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            printf("%.10g    ", D[i][k]);
        }
        printf("\n");
    }

    return D[n][n];
}

double Milne(FILE * ptr) {
    double D[n+1][n+1], h;
    int N;     

    fprintf(ptr, "Milne:\n");

    for(int i=0; i<n+1; i++)
    {
        D[i][0] = 0.;
        N = pow(2,i+2);
        h = (1.0 - 0)/(double)N;
        for(int j=0; j<=N/4 - 1; j++)
            D[i][0] += 4*h/90.*(7*f(4*j*h) + 32*f((4*j+1)*h) + 12*f((4*j+2)*h) + 32*f((4*j+3)*h) + 7*f((4*j+4)*h));
        for(int k=1; k<i+1; k++)
            D[i][k] = (pow(4,k)*D[i][k-1] - D[i-1][k-1])/(pow(4,k) - 1);
        fprintf(ptr,"%d\t%15.10g\t%15.10g\n", i, D[i][0], D[i][i]);
    }

    printf("\n\n------------------------Milne------------------------\n");
    printf("D[w][0] \t D[w][w] \n");
    for(int i =0; i <n ;i++)
    {
            printf("%.10g \t %.10g \n", D[i][0], D[i][i]);
    }
    
    printf("\n\n");

    for(int i = 0; i< n; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            printf("%.10g    ", D[i][k]);
        }
        printf("\n");
    }


    return D[n][n];
}

/*
void PrintToScreen(int size, double TAB[][size])
{
    printf("D[w][0] \t D[w][w] \n");
    for(int i =0; i <size ;i++)
    {
            printf("%.10g \t %.10g \n", TAB[i][0], TAB[i][i]);
    }
    
    printf("\n\n");

    for(int i = 0; i< size; i++)
    {
        for(int k = 0; k < i+1; k++)
        {
            printf("%.10g    ", TAB[i][k]);
        }
        printf("\n");
    }
}
*/