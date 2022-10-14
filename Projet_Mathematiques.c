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
// - Allocation d'une matrice de taille N * M
//
double** Allocation_Matrice( int N,  int M)
{
    double** Matrice = (double**)malloc(N * sizeof(double));

    for ( int i = 0; i < N; i++)
    {
        Matrice[i] = malloc(M * sizeof(double));
    }

    return Matrice;

}

//
// - Allocation d'un vecteur de taille N.
//
double* Allocation_Vecteur(int N)
{

    return  malloc(N * sizeof(double));
}

void Desallocation_Matrice(double** Matrice,  int N)
{
    for ( int i = 0; i < N; i++)
    {
        free(Matrice[i]);
    }
}

void Desallocation_Vecteur(double* Vecteur)
{
    free(Vecteur);
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                       Méthode triangulaire                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// - Algorithme de la descente
//
double* Sol_Inf(double** a, double* b,  int Taille)
{
    // - Allocation de la matrice pour X.
    double* x = Allocation_Vecteur(Taille);

    // - 
    for ( int i = 0; i < Taille; i++)
    {
        x[i] = 0;
    }

    x[0] = b[0] / a[0][0];

    for ( int i = 1; i < Taille; i++)
    {
        double Sum = 0;
        for ( int j = 0; j < i; j++)
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
double* Sol_Sup(double** a, double* b, int Taille)
{
    double* x = Allocation_Vecteur(Taille);

    for (int i = 0; i < Taille; i++)
    {
        x[i] = 0;
    }

    x[Taille - 1] = b[Taille - 1] / a[Taille - 1][Taille - 1];

    for (int i = Taille - 2; i >= 0; i--)
    {
        // - Somme
        double Sum = 0;
        for (int j = i + 1; j < Taille; j++)
        {
            Sum += a[i][j] * x[j];
        }

        x[i] = (b[i] - Sum) / a[i][i];
    }

    return x;
}

double* Gauss(double** A, double* B, int Taille)
{
    printf("Appel Gauss!!!!!!!!!!!!!!");
    for ( int i = 1; i < Taille; i++)
    {
        printf("i : %u", i);
        for ( int k = i + 1; i <= Taille; i++)
        {
            double C = A[k][i] / A[i][i];
            for ( int j = 1; i <= Taille; i++)
            {
                A[k][i] = A[k][j] - (C * A[i][j]);
            }
            B[k] = B[k] - (C * B[i]);
        }
    }

    return Sol_Sup(A, B, Taille);
}


int main()
{

    int Taille = 3;
    
    // Allocation des matrices A et B.
    double** A = Allocation_Matrice(Taille, Taille);
    double* B = Allocation_Vecteur(Taille);
        double* X;

    // - Sol inf

    // Remplissage des matrices et vecteurs.
    A[0][0] = 1;    A[0][1] = 0;    A[0][2] = 0;
    A[1][0] = 2;    A[1][1] = 3;    A[1][2] = 0;
    A[2][0] = 1;    A[2][1] = 4;    A[2][2] = -1;

    B[0] = 1;
    B[1] = 8;
    B[2] = 10;
    


    X = Sol_Inf(A, B, Taille);


    printf("La solution x pour la remontée est :\n%f\n %f\n%f\n", X[0], X[1], X[2]);
  

    // Remplissage des matrices et vecteurs.
    A[0][0] = 1;    A[0][1] = 2;    A[0][2] = 3;
    A[1][0] = 0;    A[1][1] = 4;    A[1][2] = 8;
    A[2][0] = 0;    A[2][1] = 0;    A[2][2] = 5;

    B[0] = 6;
    B[1] = 16;
    B[2] = 15;

    X = Sol_Sup(A, B, Taille);

    printf("La solution x pour la descente est :\n%f\n%f\n%f\n", X[0], X[1], X[2]);


    // - Gauss
    
    // Remplissage des
    A[0][0] = 3;    A[0][1] = 1;    A[0][2] = 2;
    A[1][0] = 3;    A[1][1] = 2;    A[1][2] = 6;
    A[2][0] = 6;    A[2][1] = 1;    A[2][2] = -1;

    
    B[0] = 2;
    B[1] = 1;
    B[2] = 4;

    X = Gauss(A, B, Taille);

    printf("La solution x pour Gauss est :\n%f\n%f\n%f\n", X[0], X[1], X[2]);


    Desallocation_Matrice(A, Taille);
    Desallocation_Vecteur(B);
    Desallocation_Vecteur(X);

    return 0;
}
