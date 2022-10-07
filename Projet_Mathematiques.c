#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          Fonctions générales                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// - Alocation d'une matrice de taille N * M
//
double** Allocation_Matrice(unsigned int N, unsigned int M)
{
    double** Matrice = (double**)malloc(N * sizeof(double));

    for (unsigned int i = 0; i < N; i++)
    {
        Matrice[i] = malloc(M * sizeof(double));
    }

    return Matrice;

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                       Méthode triangulaire                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// - Algorythme de la descente
//
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

//
// - Algorythme de la remontée
//
double* Sol_Sup(double** a, double* b, unsigned int Taille)
{
    double* x = malloc(Taille * sizeof(double));

    for (unsigned int i = 0; i < Taille; i++)
    {
        x[i] = 0;
    }

    x[Taille - 1] = b[Taille - 1] / a[Taille - 1][Taille - 1];

    for (int i = Taille - 2; i >= 0; i--)
    {
        printf("i = %u\n", i);
        double Sum = 0;
        for (unsigned int j = i + 1; j <= Taille; j++)
        {
            printf("j = %u \n", j);
            Sum = Sum + a[i][j] * x[j];
        }
        x[i] = (b[i] - Sum) / a[i][i];
    }

    return x;

}

int main()
{

    unsigned int Taille = 3;
    
    // Allocation des matrices A et B.
    double** A = Allocation_Matrice(Taille, Taille);
    double* B = (double*)malloc(Taille * sizeof(double));

    // Remplissage des
    A[0][0] = 1;
    A[0][1] = 0;
    A[0][2] = 0;
    A[1][0] = 2;
    A[1][1] = 3;
    A[1][2] = 0;
    A[2][0] = 1;
    A[2][1] = 4;
    A[2][2] = -1;

    B[0] = 1;
    B[1] = 8;
    B[2] = 10;
    
    double* X;

    X = Sol_Inf(A, B, Taille);

    printf("La solution x pour la remontée est :\n%f\n %f\n%f\n", X[0], X[1], X[2]);
  
    // Remplissage des
    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = 3;
    A[1][0] = 0;
    A[1][1] = 4;
    A[1][2] = 8;
    A[2][0] = 0;
    A[2][1] = 0;
    A[2][2] = 5;

    B[0] = 6;
    B[1] = 16;
    B[2] = 15;

    X = Sol_Sup(A, B, Taille);

    printf("La solution x pour la descente est :\n%f\n%f\n%f\n", X[0], X[1], X[2]);

    return 0;

}
