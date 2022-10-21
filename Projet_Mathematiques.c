#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          Fonctions générales                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void Nettoyer_Matrice(double** A, int N, int M)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            A[i][j] = 0;
        }
    }
}

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

void Afficher_Matrice(double**  Matrice, int Taille)
{
    for (int i = 0; i < Taille; i++)
    {
        printf("\n");

        for (int j = 0; j < Taille; j++)
        {
            printf("| %f ",Matrice[i][j]);
        }
        printf("|\n");
    }

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

void Gauss(double** A, double* B, double** U, double* e, int Taille)
{
    for (int i = 0; i < Taille - 1; i++)
    {
       for ( int k = i + 1; k < Taille; k++)
        {
            double C = A[k][i] / A[i][i];
            for ( int j = 0; j < Taille; j++)
            {
                U[k][j] = A[k][j] - (C * A[i][j]);
            }
           e[k] = B[k] - (C * B[i]);
        }
    }
}

double* Resol_Gauss(double** A, double* B, int Taille)
{
    Gauss(A,B,A,B,Taille);

    return Sol_Sup(A,B,Taille);
}

void LU(double** L, double** U, double** A, int Taille)
{
    Nettoyer_Matrice(L, Taille, Taille);
    Nettoyer_Matrice(U, Taille, Taille);

    // Rempli les 1 en diagonale
    for (int i = 0; i < Taille; i++)
    {
        L[i][i] = 1;
    }
    
    // Calcule L et U
    for (int i = 0; i < Taille - 1; i++)
    {
       for ( int k = i + 1; k < Taille; k++)
        {
            double C = A[k][i] / A[i][i];

            L[k][i] = C;
            
            for (int j = 0; j < Taille; j++)
            {
                U[k][j] = A[k][j] - (C * A[i][j]);
            }
        }
    }

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

    X = Resol_Gauss(A, B, Taille);
    
    printf("\nLa solution x pour Gauss est :\n%f\n%f\n%f\n", X[0], X[1], X[2]);
    
    // - LU
    A[0][0] = 3;    A[0][1] = 1;    A[0][2] = 2;
    A[1][0] = 3;    A[1][1] = 2;    A[1][2] = 6;
    A[2][0] = 6;    A[2][1] = 1;    A[2][2] = -1;

    double** L = Allocation_Matrice(Taille, Taille);
    double** U = Allocation_Matrice(Taille, Taille);

    LU(L, U, A, Taille);

    printf("U : \n");
    Afficher_Matrice(U, Taille);

    printf("L : \n");
    Afficher_Matrice(L, Taille);

    Desallocation_Matrice(A, Taille);
    Desallocation_Vecteur(B);
    Desallocation_Vecteur(X);

    return 0;
}
