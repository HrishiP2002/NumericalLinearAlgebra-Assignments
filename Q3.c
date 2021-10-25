#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double** Augmented[9];



double** Matrix_create(int n)
{
    double** Matrix = (double**) malloc (n * sizeof(double*));

    for(int i=0; i<n; i++)
    {
        Matrix[i] = (double*) malloc (n * sizeof(double));
    }

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            Matrix[i][j] = 1/(float)(i + j + 1); 
        }
    }

    return Matrix;
}

void gaussElimination(int m, int n, double** a){
    int i,j,k;
    for(i=0;i<m-1;i++)
    {
        //Partial Pivoting
        for(k=i+1; k<m; k++)
        {
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if(fabs(a[i][i])<fabs(a[k][i]))
            {
                //Swap the rows
                for(j=0;j<n;j++)
                {      
                    double temp;
                    temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
            }
        }
        //Begin Gauss Elimination
        for(k=i+1;k<m;k++)
        {
            double  term=a[k][i]/ a[i][i];
            for(j=0;j<n;j++)
            {
                a[k][j]=a[k][j]-term*a[i][j];
            }
        }
    }
             
}

void Augmenter(double** A, double* b, int n)
{
    double** H = (double**) malloc(n * sizeof(double*));

    for(int i=0; i<n; i++)
    {
        H[i] = (double*) malloc((n+1) * sizeof(double));
    }

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n+1; j++)
        {
            if(j != n)
            {
                H[i][j] = A[i][j];
            }
            else
            {
                H[i][j] = b[i];
            }
        }
    }


    // for(int i=0; i<n; i++)
    // {
    //     for(int j=0; j<n+1; j++)
    //     {
    //          printf("%lf\n", H[i][j]); 

    //     }
    // }

    Augmented[n-2] = H;

}

double* BackwardSubstitution(double** U , double* y, int n)  // Ux = y
{
    double* x = (double*) malloc( n* sizeof(double));

    for(int i = n-1; i>=0; i--)
    {

        double sum = 0;

        for(int j=i+1; j<n; j++)
        {
            sum = sum + U[i][j] * x[j];
        }

        x[i] = (y[i] - sum)/U[i][i];

    }

    return x;
 

}

double InfinityNorm(double* x, int n)
{
    double* y = (double*) malloc(n * sizeof(double));
    double Max = -1;

    for(int i=0; i<n; i++)
    {
        y[i] = fabs(x[i]-1);
    }

    for(int i=0; i<n; i++)
    {
        if(y[i]>Max)
        {
            Max = y[i];
        }
    }

    return Max;
}


int main()
{
    double*** Matrices = (double***) malloc (9* sizeof(double**));

    for(int i=0; i<9; i++)
    {
        Matrices[i] = Matrix_create(i+2);
    }

    double** b = (double**) malloc (9 * sizeof(double*));

    for(int i=0; i<9; i++)
    {
        b[i] = (double*) malloc (9 * sizeof(double));

        for(int j=0; j<i+2; j++)
        {
            b[i][j] = 0;

            for(int k=1; k<=i+2; k++)
            {
                b[i][j] += 1/(float)(j + k);
            }

        }
    }


    // for(int i=0; i<3; i++)
    // {
    //     for(int j=0; j<3; j++)
    //     {
    //         printf("%lf\n", Matrices[1][i][j]);
    //     }
    // }




    for(int k=0; k<9; k++)
    {

        Augmenter(Matrices[k], b[k], k+2);

        gaussElimination(k+2, k+3, Augmented[k]);
    }


    for(int h=0; h<8; h++)
    {
        for(int i=0; i<h+2; i++)
        {
            for(int j=0; j<h+3; j++)
            {
                if(j != h+3)
                {
                    Matrices[h][i][j] = Augmented[h][i][j];
                }
                else
                {
                    b[h][i] = Augmented[h][i][j]; 
                }
            }
        }
    }


    double** x = (double**) malloc(9 * sizeof(double*));

    double Errors[9];

    printf("\nThe Errors in the computed solution are as follows: \n");

    for(int i=0; i<9; i++)
    {
        x[i] = (double*) malloc((i+2) * sizeof(double));

        x[i] = BackwardSubstitution(Matrices[i], b[i], i+2);

        Errors[i] = InfinityNorm(x[i], i+2);
        printf("\nThe error for n = %d is: %lf", i+2, Errors[i]);
    }

    return 0;

    


}