#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          Fonctions générales                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// - Gestion des vecteurs et matrices

// - - Nettoyage

// Fonction qui nettoie un vecteur (remplissage par des 0).
void Nettoyage_Vecteur(double *Vecteur, int Taille)
{
    // Itère parmis les éléments du vecteur
    for (int i = 0; i < Taille; i++)
    {
        // Remplissage par des 0
        Vecteur[i] = 0;
    }
}

// - Fonction qui nettoie une matrice (remplissage par des 0).
void Nettoyage_Matrice(double **A, int Taille)
{
    // Itère parmis la première dimension
    for (int i = 0; i < Taille; i++)
    {
        // Nettoyage de la deuxième dimension
        Nettoyage_Vecteur(A[i], Taille);
    }
}

// - - Allocation

// Fonction qui alloue un vecteur de taille N.
double *Allocation_Vecteur(int Taille)
{
    // Allocation d'un espace de mémoire de taille = Taille Vecteur * Taille du type (8 octets).
    double *Vecteur = (double *)malloc(Taille * sizeof(double));

    // Nettoyage du vecteur.
    Nettoyage_Vecteur(Vecteur, Taille);

    return Vecteur;
}

// Fonction qui alloue une matrice carrée de taille N * M
double **Allocation_Matrice(int Taille)
{
    // Allocation de la première dimension
    double **Matrice = (double **)malloc(Taille * sizeof(double));

    // Allocation de la deuxième dimension
    for (int i = 0; i < Taille; i++)
    {
        Matrice[i] = Allocation_Vecteur(Taille);
    }

    return Matrice;
}

// - - Désallocation

// Fonction qui désalloue un vecteur de taille N.
void Desallocation_Vecteur(double *Vecteur)
{
    free(Vecteur);
}

// Fonction qui désalloue une matrice carrée de taille N.
void Desallocation_Matrice(double **Matrice, int Taille)
{
    for (int i = 0; i < Taille; i++)
    {
        Desallocation_Vecteur(Matrice[i]);
    }
}

// - - Affichage

// Fonction qui affiche un vecteur de taille N.
void Afficher_Vecteur(double *Vecteur, int Taille)
{
    // Itère parmis les éléments du vecteur.
    for (int i = 0; i < Taille; i++)
    {
        // Affichage de la valeur suivi d'un saut de ligne.
        printf("%f\n", Vecteur[i]);
    }
}

// Fonction qui affiche une matrice carrée de taille N.
void Afficher_Matrice(double **Matrice, int Taille)
{
    // Itère parmis la première dimension de la matrice.
    for (int i = 0; i < Taille; i++)
    {

        // Itère parmis la deuxième dimension de la matrice.
        for (int j = 0; j < Taille; j++)
        {
            // Affichage de la valeur avec un séparateur.
            printf("| %f ", Matrice[i][j]);
        }
        // Retour à la ligne.
        printf("|\n");
    }
}

// - Méthodes directes de résolution.

// - - Méthode triangulaire

// Fonction qui résoud une matrice carrée triangulaire inférieure (algorithme de la redescente).
double *Sol_Inf(double **A, double *B, int Taille)
{
    // Allocation du vecteur X.
    double *X = Allocation_Vecteur(Taille);

    // Calcul du premier terme de X.
    X[0] = B[0] / A[0][0];

    // Itère parmis les lignes de la matrice.
    for (int i = 1; i < Taille; i++)
    {
        // Calcul de la somme des a[i][j] * x[j]
        double Sum = 0;
        for (int j = 0; j < i; j++)
        {
            Sum = Sum + A[i][j] * X[j];
        }
        // Calcul du terme X[i].
        X[i] = (B[i] - Sum) / A[i][i];
    }

    return X;
}

// Fonction qui résoud une matrice carrée triangulaire supérieure (algorithme de la remontée).
double *Sol_Sup(double **A, double *B, int Taille)
{
    // Allocation du vecteur X.
    double *X = Allocation_Vecteur(Taille);

    // Calcul du dernier terme de X (terme initial).
    X[Taille - 1] = B[Taille - 1] / A[Taille - 1][Taille - 1];

    // Itère parmis les lignes de la matrice.
    for (int i = Taille - 2; i >= 0; i--)
    {
        // Calcul de la somme des a[i][j] * x[j]
        double Sum = 0;
        for (int j = i + 1; j < Taille; j++)
        {
            Sum += A[i][j] * X[j];
        }

        X[i] = (B[i] - Sum) / A[i][i];
    }

    return X;
}

// - - Méthode de Gauss

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
    Nettoyage_Matrice(L, Taille);

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

// - Méthodes itératives de résolution.

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

    // Remplissage de X_k par des 1.
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

    // Les différentes parties sont mises entre crochet afin que les variables aients une portée locale (propres à chaques parties).

    // - Méthodes directes

    // - - Méthode triangulaire inférieure.

    {
        // Définition de la taille.
        int Taille = 3;
        // Allocation des matrices et vecteurs.
        double **A = Allocation_Matrice(Taille);
        double *B = Allocation_Vecteur(Taille);

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

        double *X = Sol_Inf(A, B, Taille);

        printf("La solution X de l'équation AX = B avec la méthode triangulaire inférieure est :\n");
        Afficher_Vecteur(X, Taille);

        // Désallocation des matrices et vecteurs.
        Desallocation_Matrice(A, Taille);
        Desallocation_Vecteur(B);
        Desallocation_Vecteur(X);
    }

    // - - Méthode triangulaire supérieure.

    {
        // Définition de la taille.
        int Taille = 3;
        // Allocation des matrices et vecteurs.
        double **A = Allocation_Matrice(Taille);
        double *B = Allocation_Vecteur(Taille);

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

        double *X = Sol_Sup(A, B, Taille);

        printf("La solution X de l'équation AX = B avec la méthode triangulaire supérieure est :\n");
        Afficher_Vecteur(X, Taille);

        // Désallocation des matrices et vecteurs.
        Desallocation_Matrice(A, Taille);
        Desallocation_Vecteur(B);
        Desallocation_Vecteur(X);
    }

    // - - Méthode de Gauss.

    {
        // Définition de la taille.
        int Taille = 3;
        // Allocation des matrices et vecteurs.
        double **A = Allocation_Matrice(Taille);
        double *B = Allocation_Vecteur(Taille);

        // Remplissage de A
        A[0][0] = 1;
        A[0][1] = 2;
        A[0][2] = 3;
        A[1][0] = 5;
        A[1][1] = 2;
        A[1][2] = 1;
        A[2][0] = 3;
        A[2][1] = -1;
        A[2][2] = 1;

        // Remplissage de B
        B[0] = 5;
        B[1] = 5;
        B[2] = 6;

        // Résolution de l'équation.
        double *X = Resol_Gauss(A, B, Taille);

        // Affichage
        printf("La solution x pour Gauss est :\n");
        Afficher_Vecteur(X, Taille);

        // Désallocation des matrices et vecteurs.
        Desallocation_Matrice(A, Taille);
        Desallocation_Vecteur(B);
        Desallocation_Vecteur(X);
    }

    // - - Résolution par factorisation LU

    {
        // Définition de la taille.
        int Taille = 3;
        // Allocation des matrices et vecteurs.
        double **A = Allocation_Matrice(Taille);
        double *B = Allocation_Vecteur(Taille);

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

        // Factorisation de A = LU
        LU(L, A, Taille);

        // Résolution de LY = B
        double *Y = Sol_Inf(L, B, Taille);

        // Résolution de UX = Y
        double *X = Sol_Sup(A, Y, Taille);

        // Affichage
        printf("La solution X pour LU : \n");
        Afficher_Vecteur(X, Taille);

        // Désallocation des matrices et vecteurs.
        Desallocation_Matrice(L, Taille);
        Desallocation_Matrice(A, Taille);
        Desallocation_Vecteur(B);
        Desallocation_Vecteur(Y);
        Desallocation_Vecteur(X);
    }

    // - - Cholesky

    {
        // Définition de la taille.
        int Taille = 4;
        // Allocation des matrices et vecteurs.
        double **A = Allocation_Matrice(Taille);
        double *B = Allocation_Vecteur(Taille);

        // Remplissage de A
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

        // Remplissage de B
        B[0] = 5;
        B[1] = 1;
        B[2] = 3;
        B[3] = 1;

        printf("La solution pour Choleski est :\n");

        // Factorisation de A = LL^T
        double** L = Cholesky(A, Taille);

        // Résolution de LY = B
        double *Y = Sol_Inf(L, B, Taille);

        // Résolution de L^TX = Y
        Transposer(L, Taille);

        double *X = Sol_Sup(L, Y, Taille);

        printf("La solution X pour Cholesky : \n");

        Afficher_Vecteur(X, Taille);

        Desallocation_Matrice(A, Taille);
        Desallocation_Matrice(L, Taille);
        Desallocation_Vecteur(Y);
        Desallocation_Vecteur(B);
        Desallocation_Vecteur(X);
    }

    // - Méthodes itératives

    // - - Méthode de Jacobi
    {   
        // Définition de la taille
        int Taille = 3;
        // Allocation des matrices et vecteurs.
        double** A = Allocation_Matrice(Taille);
        double* B = Allocation_Vecteur(Taille);

        // Remplissage de A
        A[0][0] = 4;
        A[0][1] = 1;
        A[0][2] = 1;
        A[1][0] = 1;
        A[1][1] = 3;
        A[1][2] = 1;
        A[2][0] = 2;
        A[2][1] = 0;
        A[2][2] = 5;
        
        // Remplissage de B
        B[0] = 1;
        B[1] = 1;
        B[2] = 1;

        double* X = Jacobi(A, B, Taille);

        printf("La solution X de AX = B avec la méthode Jacobi : \n");

        Afficher_Vecteur(X, Taille);

        // - Desallocation avant la fermeture du programme.

        Desallocation_Vecteur(B);
        Desallocation_Vecteur(X);
        Desallocation_Matrice(A, Taille);
    }

    return 0;
}
