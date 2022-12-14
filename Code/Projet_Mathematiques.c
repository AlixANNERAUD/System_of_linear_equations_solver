#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          Fonctions générales                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void Nettoyer_Matrice(double **A, int Taille)
{
    // Itère parmis les premières dimension
    for (int i = 0; i < Taille; i++)
    {
        // Itère parmis les deuxième dimension
        for (int j = 0; j < Taille; j++)
        {
            A[i][j] = 0;
        }
    }
}

//
// - Allocation d'un vecteur de taille N.
//
double *Allocation_Vecteur(int Taille)
{
    // - Allocation d'un espace de mémoire de taille = Taille Vecteur * Taille du type (8 octets).
    return malloc(Taille * sizeof(double));
}

//
// - Allocation d'une matrice carrée de taille N * M
//
double **Allocation_Matrice(int Taille)
{
    // - Allocation de la première dimension
    double **Matrice = (double **)malloc(Taille * sizeof(double));

    // - ALlocation de la deuxième dimension
    for (int i = 0; i < Taille; i++)
    {
        Matrice[i] = Allocation_Vecteur(Taille);
    }
    
    return Matrice;
}



void Desallocation_Matrice(double **Matrice, int Taille)
{
    for (int i = 0; i < Taille; i++)
    {
        free(Matrice[i]);
    }
}

void Desallocation_Vecteur(double *Vecteur)
{
    free(Vecteur);
}

void Afficher_Matrice(double **Matrice, int Taille)
{
    for (int i = 0; i < Taille; i++)
    {
        printf("\n");

        for (int j = 0; j < Taille; j++)
        {
            printf("| %f ", Matrice[i][j]);
        }
        printf("|\n");
    }
}

void Afficher_Vecteur(double *Vecteur, int Taille)
{
    for (int i = 0; i < Taille; i++)
    {
        printf("%f\n", Vecteur[i]);
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
double *Sol_Inf(double **a, double *b, int Taille)
{
    // - Allocation de la matrice pour X.
    double *x = Allocation_Vecteur(Taille);

    // -
    for (int i = 0; i < Taille; i++)
    {
        x[i] = 0;
    }

    x[0] = b[0] / a[0][0];

    for (int i = 1; i < Taille; i++)
    {
        double Sum = 0;
        for (int j = 0; j < i; j++)
        {
            Sum = Sum + a[i][j] * x[j];
        }

        x[i] = (b[i] - Sum) / a[i][i];
    }

    return x;
}
//
// - Algorithme de la remontée
//
double *Sol_Sup(double **a, double *b, int Taille)
{
    double *x = Allocation_Vecteur(Taille);

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

void Gauss(double **A, double *B, double **U, double *e, int Taille)
{
    for (int i = 0; i < Taille - 1; i++)
    {
        for (int k = i + 1; k < Taille; k++)
        {
            double C = A[k][i] / A[i][i];
            for (int j = 0; j < Taille; j++)
            {
                U[k][j] = A[k][j] - (C * A[i][j]);
            }
            e[k] = B[k] - (C * B[i]);
        }
    }
}

double *Resol_Gauss(double **A, double *B, int Taille)
{
    Gauss(A, B, A, B, Taille);

    return Sol_Sup(A, B, Taille);
}

// Transforme A en U.
void LU(double **L, double **A, int Taille)
{
    Nettoyer_Matrice(L, Taille);

    // Rempli les 1 en diagonale
    for (int i = 0; i < Taille; i++)
    {
        L[i][i] = 1;
    }

    // Calcule L et U
    for (int i = 0; i < Taille - 1; i++)
    {
        for (int k = i + 1; k < Taille; k++)
        {
            double C = A[k][i] / A[i][i];

            L[k][i] = C;

            for (int j = 0; j < Taille; j++)
            {
                A[k][j] = A[k][j] - (C * A[i][j]);
            }
        }
    }
}

// Fonction qui transforme A en L*L^T. Renvoi L.
double **Cholesky(double **A, int Taille)
{
    double **L = Allocation_Matrice(Taille); // Allocation de la matrice L.
    Nettoyer_Matrice(L, Taille);

    for (int j = 0; j < Taille; j++)
    {
        double Sum = 0;
        // Calcul de la diagonale
        // Calcul de la somme

        for (int k = 0; k < j; k++)
        {

            Sum += L[j][k] * L[j][k];
        }

        // Calcul des termes de la diagonale de L.
        L[j][j] = sqrt(A[j][j] - Sum);

        // Calcul des termes pour i = j+1 à N.
        for (int i = j + 1; i < Taille; i++)
        {

            Sum = 0;
            for (int k = 0; k < j; k++)
            {

                Sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - Sum) / L[j][j];
        }
    }
    return L;
}

// Fonction qui transpose une matrice carrée inférieure en matrice carrée supérieure.
double **Transposer(double **Matrice, int Taille)
{
    // On itère parmis les colones.
    for (int i = 0; i < Taille; i++)
    {
        // On itère parmis les lignes.
        for (int j = i + 1; j < Taille; j++)
        {
            // On inverse les termes (en cooordonnées) de la matrice.
            Matrice[i][j] = Matrice[j][i];
            // On supprime les termes inférieures.
            Matrice[j][i] = 0;
        }
    }
    return Matrice;
}

double *Jacobi(double **A, double *B, int Taille)
{
    double *X_k = Allocation_Vecteur(Taille);
    for (int i = 0; i < Taille; i++)
    {
        X_k[i] = 1;
    }
    double *X_k_1 = Allocation_Vecteur(Taille);

    double Norme = 1;

    int it = 0;
    int it_max = 10;

    double Epsilon = 0.001;

    while ((Norme > Epsilon) && (it < it_max))
    {
        // Calcul de X_k_1
        for (int i = 0; i < Taille; i++)
        {
            // Calcul de la somme
            double Somme = 0;
            for (int j = 0; j < Taille; j++)
            {
                if (i != j)
                {
                    Somme += A[i][j] * X_k[j];
                }
            }
            X_k_1[i] = (B[i] - Somme) / A[i][i];
        }
        // Calcul de la norme
        Norme = 0;
        for (int i = 0; i < Taille; i++)
        {
            Norme += (X_k_1[i] - X_k[i]) * (X_k_1[i] - X_k[i]);
        }
        Norme = sqrt(Norme);

        for (int i = 0; i < Taille; i++)
        {
            X_k[i] = X_k_1[i];
        }
        it++;
    }
    return X_k_1;
}

int main()
{
    int Taille = 3;

    // Allocation des matrices A et B.
    double **A = Allocation_Matrice(Taille);
    double *B = Allocation_Vecteur(Taille);
    double *X;

    // - Sol inf

    // Remplissage des matrices et vecteurs.
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

    X = Sol_Inf(A, B, Taille);

    printf("La solution x pour la remontée est :\n%f\n %f\n%f\n", X[0], X[1], X[2]);

    // Remplissage des matrices et vecteurs.
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

    // - Gauss

    // Remplissage de la matrice
    A[0][0] = 3;
    A[0][1] = 1;
    A[0][2] = 2;
    A[1][0] = 3;
    A[1][1] = 2;
    A[1][2] = 6;
    A[2][0] = 6;
    A[2][1] = 1;
    A[2][2] = -1;

    B[0] = 2;
    B[1] = 1;
    B[2] = 4;

    X = Resol_Gauss(A, B, Taille);

    printf("\nLa solution x pour Gauss est :\n%f\n%f\n%f\n", X[0], X[1], X[2]);

    // - LU
    Nettoyer_Matrice(A, Taille);

    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = 3;
    A[1][0] = 5;
    A[1][1] = 2;
    A[1][2] = 1;
    A[2][0] = 3;
    A[2][1] = -1;
    A[2][2] = 1;

    B[0] = 5;
    B[1] = 5;
    B[2] = 6;

    double **L = Allocation_Matrice(Taille);

    LU(L, A, Taille);

    double *Y = Sol_Inf(L, B, Taille);

    X = Sol_Sup(A, Y, Taille);

    printf("La solution X pour LU : \n");
    Afficher_Vecteur(X, Taille);

    Desallocation_Matrice(L, Taille);

    Desallocation_Matrice(A, Taille);

    Desallocation_Vecteur(Y);

    Desallocation_Vecteur(X);

    Desallocation_Vecteur(B);

    // - Cholesky

    Taille = 4;

    A = Allocation_Matrice(Taille);
    B = Allocation_Vecteur(Taille);

    Y = Allocation_Vecteur(Taille);
    X = Allocation_Vecteur(Taille);

    Nettoyer_Matrice(A, Taille);

    A[0][0] = 1;
    A[0][1] = 1;
    A[0][2] = 1;
    A[0][3] = 1;
    A[1][0] = 1;
    A[1][1] = 5;
    A[1][2] = 5;
    A[1][3] = 5;
    A[2][0] = 1;
    A[2][1] = 5;
    A[2][2] = 14;
    A[2][3] = 14;
    A[3][0] = 1;
    A[3][1] = 5;
    A[3][2] = 14;
    A[3][3] = 15;

    B[0] = 5;
    B[1] = 1;
    B[2] = 3;
    B[3] = 1;

    Afficher_Matrice(A, Taille);

    printf("La solution pour Choleski est :\n");

    L = Cholesky(A, Taille);

    printf("La matrice L : \n");

    Afficher_Matrice(L, Taille);

    Y = Sol_Inf(L, B, Taille);

    Transposer(L, Taille);

    printf("La matrice transaposée est : \n");

    Afficher_Matrice(L, Taille);

    X = Sol_Sup(L, Y, Taille);

    printf("La solution X pour Cholesky : \n");

    Afficher_Vecteur(X, Taille);

    Desallocation_Matrice(L, Taille);

    Desallocation_Vecteur(Y);

    Desallocation_Matrice(A, Taille);

    Desallocation_Vecteur(B);

    Desallocation_Vecteur(X);

    // - Jacobi

    Taille = 3;

    A = Allocation_Matrice(Taille);
    B = Allocation_Vecteur(Taille);

    A[0][0] = 4;
    A[0][1] = 1;
    A[0][2] = 1;
    A[1][0] = 1;
    A[1][1] = 3;
    A[1][2] = 1;
    A[2][0] = 2;
    A[2][1] = 0;
    A[2][2] = 5;

    B[0] = 1;
    B[1] = 1;
    B[2] = 1;

    X = Jacobi(A, B, Taille);

    printf("La solution X pour Jacobi : \n");

    Afficher_Vecteur(X, Taille);

    // - Desallocation avant la fermeture du programme.

    Desallocation_Vecteur(B);
    Desallocation_Vecteur(X);
    Desallocation_Matrice(A, Taille);

    return 0;
}
