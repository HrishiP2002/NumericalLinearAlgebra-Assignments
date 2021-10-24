#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void gaussElimination(int m, int n, double a[5][6]){
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

void BackwardSubstitution(double U[5][5], double  y[5], int n)  // Ux = y
{
    double x[5];

    for(int i = n-1; i>=0; i--)
    {

        double sum = 0;

        for(int j=i+1; j<n; j++)
        {
            sum = sum + U[i][j] * x[j];
        }

        x[i] = (y[i] - sum)/U[i][i];



    }
 
         printf("\n\nThe solution of the given system of Linear equations is:\n");

        for(int i=0; i<n; i++)
        {
            printf("%lf\n", x[i]);
        }

}


int main()
{
    double A[5][5];
    double x[5] , b[5];

    printf("Please enter the elements of the matrix:\n ");

    for(int i=0; i<5; i++)
    {
        for(int j=0; j<5; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("\n\nPlease enter the elements of the RHS vector: \n");

    for(int i=0; i<5; i++)
    {
        scanf("%lf", &b[i]);
    }

    double H[5][6];       // H is the augmented matrix

    for(int i=0; i<5; i++)  // copying matrix A and b into H
    {
        for(int j=0; j<5; j++)
        {
            H[i][j] = A[i][j];
        }
    }

    for(int k=0; k<5; k++)
    {
        H[k][5] = b[k];
    }

    gaussElimination(5, 6, H);   // Gauss elimination of H

    
    for(int i=0; i<5; i++)   // Copying back the right side of the augmented matrix into b
    {
        b[i] = H[i][5];
    }

    for(int i=0; i<5; i++)  // Copying back the left side of the augmented matrix into A
    {
        for(int j=0; j<5; j++)
        {
            A[i][j] = H[i][j];
        }
    }

    
    BackwardSubstitution(A, b, 5);




    
}

//Matrix A

// 2 
// 10 
// 8 
// 8 
// 6
// 1 
// 4 
// −2 
// 4 
// −1
// 0 
// 2 
// 3 
// 2 
// 1
// 3 
// 8 
// 3 
// 10 
// 9
// 1 
// 4 
// 1 
// 2 
// 1


//Vector b

// 52
// 14
// 12
// 51
// 15


//Vector c

// 50
// 4
// 12
// 48
// 12
