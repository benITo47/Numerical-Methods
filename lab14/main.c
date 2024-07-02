#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/gammp.c"
#include "/home/NR/numerical_recipes.c/gcf.c"
#include "/home/NR/numerical_recipes.c/gammln.c"
#include "/home/NR/numerical_recipes.c/gser.c"
#include "/home/NR/numerical_recipes.c/erff.c"


#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define N 10000
#define K 12
double*  generateUniform(FILE* ptr,  double* tab, long a, int c, long m, FILE* ptr2);
void calculateHistogram(double* x, FILE* ptr, double xmin, double xmax);
double normalFunc(double x, double u, double sigma);
double* generateNormal(FILE* ptr, double* x);

int main(void) {
    FILE* fptr1 = fopen("U.dat", "w");
    FILE* fptr2 = fopen("U_hist.dat", "w");
    FILE* fptr3 = fopen("N_hist.dat", "w");
    
    long a1 = 123;
    long a2 = 69069;

    int c1 = 1;
    int c2 = 1;
    
    long int m1 = pow(2, 15);
    long int m2 = pow(2, 32);

     double tab1[N];
     double tab2[N];
    double normal_tab[N];

    double* tabl1 = generateUniform(fptr1, tab1, a1, c1, m1, fptr2);
    double* tabl2 = generateUniform(fptr1, tab2, a2, c2, m2, fptr2);
    double* tabl3 = generateNormal(fptr3,tabl2);

    

    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);

    return 0;
}

double* generateUniform(FILE* ptr,  double* x, long a, int c, long m, FILE* ptr2) {
    long X[N];
    X[0] = 10.0; 
    long double srednia = 0.0; 
    for (int i = 1; i < N; i++) {
        X[i] = (a * X[i-1] + c) % m;
        x[i] = X[i]/(m+1.0);
        srednia += x[i];
    }
    
    for (int i = 0; i < N-1; i++) {
        fprintf(ptr, "%g %g\n", x[i], x[i+1]);
    }
    fprintf(ptr, "\n\n\n\n");

    srednia /= N; 
    double sigma = 0.0; 
    for(int i = 0 ; i < N ; i++)
    {
        sigma += (x[i] - srednia)*(x[i] - srednia);
    }
    sigma = sigma/N;
    sigma = sqrt(sigma);
    //Sigma i srednia sa tutaj juz policzone 
    calculateHistogram(x,ptr2,0,1);
    return x; 

}


void calculateHistogram(double* x, FILE* ptr, double xmin, double xmax) {
    int bins[K] = {0};
    double delta = (xmax - xmin) / K;

    for (int i = 0; i < N; i++) {
        int j = (x[i] - xmin) / delta;
        if (j >= 0 && j < K) {
            bins[j]++;
        }
    }

    for (int j = 0; j < K; j++) {
        double bin_center = xmin + (j + 0.5) * delta;
        fprintf(ptr, "%g %g\n", bin_center, (double)bins[j] / N);
    }
    fprintf(ptr, "\n\n\n\n\n\n");
}

double* generateNormal(FILE* ptr, double* u2) {
    double avg = 0.2; 
    double sigma = 0.5; 
    double xmin = avg - 3 * sigma; 
    double xmax = avg + 3 * sigma; 
    double x[N];
    int count = 0;

    while (count < N) {
        // Ensure we have pairs of uniform random numbers from u2 array
        for (int i = 0; i < N - 1; i += 2) {
            double u1 = u2[i];
            double u2_val = u2[i + 1];

            // Scale u1 to the desired range
            double y = xmin + (xmax - xmin) * u1;
            double d = normalFunc(y, avg, sigma);

            // Perform rejection sampling
            if (u2_val <= d) {
                x[count] = y;
                count++;
                if (count >= N) {
                    break;
                }
            }
        }
    }


    double srednia = 0.0; 
    for (int i = 0; i < N; i++) {
        srednia += x[i];
    }
    srednia = srednia / N;


    double sigmaCalc = 0.0; 
    for (int i = 0; i < N; i++) {
        sigmaCalc += pow((x[i] - srednia), 2); 
    }
    sigmaCalc = sigmaCalc / N; 
    sigmaCalc = sqrt(sigmaCalc);


    printf("Dla rozkladu normalnego:\n");
    printf("Srednia: %f \t Odchylenie standardowe: %f\n", srednia, sigmaCalc);


    printf("Teoretyczne srednia: %f \t Teoretyczne odchylenie standardowe: %f\n", avg, sigma);

    

    calculateHistogram(x, ptr, xmin, xmax);
    return x; 
}

double normalFunc(double x, double u, double sigma) {
    return (1 / (sigma * sqrt(2 * M_PI))) * exp(-pow((x - u), 2) / (2 * sigma * sigma));
}
