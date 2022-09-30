#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

double** Allocation_Matrice(unsigned int N, unsigned int M)
{
    double** Matrice = (double**)malloc(N * sizeof(double));

    for (unsigned int i = 0; i < N; i++)
    {
        Matrice[i] = malloc(M * sizeof(double));
    }

    return Matrice;

}

double* Sol_Inf(double** a, double* b, unsigned int Taille)
{

    double* x = malloc(Taille*sizeof(double));

    for (unsigned int i = 0; i < Taille; i++)
    {
        x[i] = 0;
    }

    x[0] = b[0] / a[0][0];

    for (unsigned int i = 1; i < Taille; i++)
    {
        double Sum = 0;
        for (unsigned int j = 0; j < i; j++)
        {
            Sum = Sum + a[i][j]*x[j];
        }
        
        x[i] = (b[i] - Sum) / a[i][i];
    }
    
    return x;
}

/*
double* algo_montee(double** a, double* b, unsigned int Taille)
{
    double* x = malloc(Taille * sizeof(double));

    for (unsigned int k = 0; k < Taille; k++)
    {
        x[Taille] = 0;
    }

}*/

int main()
{

    unsigned int Taille = 3;
    
    // Allocation des matrices.
    double** A = Allocation_Matrice(Taille, Taille);

    A[0][0] = 1;
    A[0][1] = 0;
    A[0][2] = 0;
    A[1][0] = 2;
    A[1][1] = 3;
    A[1][2] = 0;
    A[2][0] = 1;
    A[2][1] = 4;
    A[2][2] = -1;

    double* B = (double*)malloc(Taille * sizeof(double));
    
    B[0] = 1;
    B[1] = 8;
    B[2] = 10;
    
    double* X;


    // Saisie de la matrice

    X = Sol_Inf(A, B, Taille);

    printf("La solution x est :\n");

    for (unsigned int i = 0; i < Taille; i++)
    {
        printf("%f\n", X[i]);
    }


    //algo_descente(a, b, Taille);

    return 0;

}