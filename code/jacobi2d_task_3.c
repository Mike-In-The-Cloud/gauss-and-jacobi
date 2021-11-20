#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

void jacobi(int m, int n, int nthreads, double tol)
{
    int iter;
    double tstart, tstop, priv_difmax; // = 0.0001;
    register int i, j; 
    float top, bottom, left, right;         // declaring variables for boundary conditions

    top = 30;
    bottom = 60;
    left = 110;
    right = 140;

     /* define the number of threads to be used */
    omp_set_num_threads(nthreads);
    double t[m + 2][n + 2], tnew[m + 1][n + 1], diff, difmax, priv_divmax;

    printf("m = %d n = %d tol = %lf threads = %d\n", m, n, tol, nthreads);

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
        t[i][0] = left;      // left
        t[i][n + 1] = right; // right
    }
    for (j = 1; j <= n; j++)
    {
        t[0][j] = top;        // top
        t[m + 1][j] = bottom; // bottom
    }

    // main loop
    iter = 0;
    difmax = 1000000.0;
    while (difmax > tol) // changing while loop for single itteration
    {
        iter++;
#pragma omp parallel default(shared) private(i, j, diff, priv_difmax)
        {
// update temperature for next iteration
            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;
                }
            }

            // work out maximum difference between old and new temperatures
            difmax = 0.0;

            priv_difmax = 0.0;
#pragma omp for schedule(static)
            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    diff = fabs(tnew[i][j] - t[i][j]);
                    if (diff > priv_difmax)
                    {
                        priv_difmax = diff;
                    }
                    // copy new to old temperatures
                    t[i][j] = tnew[i][j];
                }
            }
#pragma omp critical
            if (priv_difmax > difmax)
            {
                difmax = priv_difmax;
            }
		}// end parallel region
    }

    tstop = omp_get_wtime();

    // print results
    printf("iter = %d  difmax = %9.11lf\n", iter, difmax);
    printf("time taken is %4.6lf\n", (tstop - tstart));
    /* for (i = 0; i <= m + 1; i++)
	{
		printf("\n");
		for (j = 0; j <= n + 1; j++)
		{
			printf("%3.5lf ", t[i][j]);
		}
	}   */
    printf("\n");
}


int main(int argc, char *argv[])
{

   
    int m, n, nthreads, i;
    double tol;
    
    
    m = atoi(argv[1]);        // itterations
    n = atoi(argv[2]);        // width of array
    tol = atof(argv[3]);      // diffmax  
    nthreads = atoi(argv[4]); // number of threads
    
    
    jacobi(m, n, nthreads, tol);
   
}