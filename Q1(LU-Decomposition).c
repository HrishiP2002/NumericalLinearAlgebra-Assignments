#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// struct LU
// {
//     int L[4][4];
//     int U[4][4];
// };

double L[4][4];
double U[4][4];

// int* innerProduct(int* vect1, int* vect2, int size)
// {
//     int* prod[size];

//     for(int i=0; i<size; i++)
//     {
//         prod[i] = vect1[i] * vect2[i];
//     }

//     return prod;
// }


void LUdecomposition(double matrix[4][4], int matrix_size)
{
    double Lo[4][4];
    double Up[4][4];

    struct LU* H;

    for(int i=0; i<matrix_size; i++)
    {
        for(int j=0; j<matrix_size; j++)
        {
            Lo[i][j] = 0;
            Up[i][j] = 0;
        }
    }



    for(int i=0; i<matrix_size; i++)
    {
        //Upper Triangular: 

        for(int k=i; k<matrix_size; k++)
        {
            double sum = 0;

            for(int j=0; j<i; j++)
                sum = sum + Lo[i][j] * Up[j][k];


            Up[i][k] = matrix[i][k] - sum;
        }

        // Lower Triangular:

        for(int k=i; k<matrix_size; k++)
        {
            if(i == k)
            {
                Lo[i][i] = 1;
            }

            else
            {
                double summ = 0;
                for(int j=0; j<i; j++)
                    summ = summ + Lo[k][j] * Up[j][i];
                
                Lo[k][i] = (matrix[k][i] - summ)/Up[i][i];
            }
        }
    }

    printf("\n\nThe unit lower triangular matrix is:\n\n ");

    for(int i=0; i<matrix_size; i++)
    {
        for(int j=0; j<matrix_size; j++)
        {
            printf("%lf\t", Lo[i][j]);
        }

        printf("\n");
    }

    printf("\n\n\n");
    printf("The upper triangular matrix is:\n\n");


    for(int i=0; i<matrix_size; i++)
    {
        for(int j=0; j<matrix_size; j++)
        {
            printf("%lf\t", Up[i][j]);
        }

        printf("\n");
    }

    // H->L = Lo;

    // H->U = Up;



    for(int i=0; i<matrix_size; i++)
    {
        for(int j=0; j<matrix_size; j++)
        {
            L[i][j] = Lo[i][j];
            U[i][j] = Up[i][j];
        }
    }



}


void BackwardSubstitution( double* y, int matrix_size)  // Ux = y
{
    double* x = (double*) malloc (matrix_size * sizeof(double));


    for(int i = matrix_size-1; i>=0; i--)
    {
        double sum = 0;

        for(int j=i+1; j<matrix_size; j++)
        {
            sum = sum + U[i][j] * x[j];
        }

        x[i] = (y[i] - sum)/U[i][i];


    }


            printf("\n\nThe solution of the given system of Linear equations is:\n");

        for(int i=0; i<matrix_size; i++)
        {
            printf("%lf\n", x[i]);
        }
}




void ForwardSubstitution( double* b, int matrix_size)  // Ly = b
{

    double* y = (double*) malloc (matrix_size * sizeof(double));

    for(int i=0; i<matrix_size; i++)
    {
        double sum = 0;

        for(int j=0; j<i; j++)
        {
            sum = sum + L[i][j] * y[j];
        }

        y[i] = b[i] - sum;
    }

    BackwardSubstitution( y , matrix_size);
}


int main()
{
    int matrix_size = 4;

    double matrix[4][4];


    double vector[4];

    printf("Please enter the entries of the matrix line by line: \n");

    for(int i=0; i<matrix_size; i++)
    {
        for(int j=0; j<matrix_size; j++)
        {
            scanf("%lf", &matrix[i][j]);
        }
    }

    printf("\n\nPlease enter the entries of the RHS vector(up to down): \n");

    for(int i=0; i<matrix_size; i++)
    {
        scanf("%lf", &vector[i]);
    }

    LUdecomposition(matrix, matrix_size);

    ForwardSubstitution( vector, matrix_size);

    return 0;

    



}


//Matrix A:

// 2
// 1
// -1
// 3
// -2
// 0
// 0
// 0
// 4
// 1
// -2
// 6
// -6
// -1
// 2
// 3


//vector b

// 13
// -2
// 24
// -14


//vector c

// 12
// -8
// 21
// 26
