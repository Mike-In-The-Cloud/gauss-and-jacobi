#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>



void gauss(int m, int n, int nthreads, double tol)
{

    int iter;
    double tstart, tstop, priv_difmax; // = 0.0001;
    register short int i, j;                // added register keyword for faster access to counters

    float top, bottom, left, right;         // declaring variables for boundary conditions

    top = 30.0;
    bottom = 60.0;
    left = 110.0;
    right = 140.0;

    /* define the number of threads to be used */
    omp_set_num_threads(nthreads);

    double t[m + 2][n + 2], tnew[m + 1][n + 1], diff, difmax;

    printf("%d %d %lf %d\n", m, n, tol, nthreads);
    tstart = omp_get_wtime();
    // initialise temperature array
    for (i = 0; i <= m + 1; i++)
    {
        for (j = 0; j <= n + 1; j++)
        {
            t[i][j] = 30.0;
        }
    }

    // fix boundary conditions
    for (i = 1; i <= m; i++)
    {
        t[i][0] = left;
        t[i][n + 1] = right;
    }
    for (j = 1; j <= n; j++)
    {
        t[0][j] = top;        // top
        t[m + 1][j] = bottom; // bottom
    }

    // main loop
    iter = 0;
    difmax = 1000000.0;

while (iter < 2) // changing while loop for single itteration
    {

        iter++;
        difmax = 0.0;
        // update temperature for next iteration
#pragma omp parallel default(shared) private(diff, priv_difmax, i, j)
        {

            priv_difmax = 0.0;
#pragma omp parallel for schedule(static)
            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;

                    // work out maximum difference between old and new temperatures
                    diff = fabs(tnew[i][j] - t[i][j]);
                    if (diff > priv_difmax)
                    {
                        priv_difmax = diff;
                    }
                    t[i][j] = tnew[i][j];
                }
            }
#pragma omp critical
            if (priv_difmax > difmax)
            {
                difmax = priv_difmax;
            }
        } // end parallel region
    }
    tstop = omp_get_wtime(); // stop omp timer
    // print results
    printf("iter = %d  difmax = %9.11lf\n", iter, difmax);

    printf("time taken is %4.6lf\n", (tstop - tstart));
    for (i = 0; i <= m + 1; i++)
    {
        printf("\n");
        for (j = 0; j <= n + 1; j++)
        {
            printf("%3.5lf ", t[i][j]);
        }
    }
    printf("\n");
}

int main(int argc, char *argv[])
{

    int m, n, nthreads, i, nloops;
    double tol;
    
    
    m = atoi(argv[1]);        // itterations
    n = atoi(argv[2]);        // width of array
    tol = atof(argv[3]);      // diffmax
    nthreads = atoi(argv[4]); // number of threads

    
    printf("Number of threads: %d\n", nthreads);
    
    gauss(m, n, nthreads, tol);
    
}