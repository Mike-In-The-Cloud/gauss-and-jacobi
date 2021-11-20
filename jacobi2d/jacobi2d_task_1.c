#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char *argv[])
{

//**********************************//
//variables for start and stop timer//
//**********************************//
struct timeval startTime, stopTime;

// total time variable
long totalTime;

struct timeVal{
    time_t	tv_sec; // seconds 
    suseconds_t tv_usec; // microseconds
};


int m; 
int n;
double tol;// = 0.0001;

float top, bottom, left, right; // declaring variables for boundary conditions

top = 30;
bottom = 60;
left = 110;
right = 140;

register short int i, j;// added register keyword for faster access to counters

int iter; 
m = atoi(argv[1]); // itterations
n = atoi(argv[2]); // width of array
tol = atof(argv[3]); // diffmax

//**********************************//
//            pointers              //
//**********************************//

float *top_pointer, *bottom_pointer, *left_pointer, *right_pointer = NULL;

// assigning addresses to pointers
top_pointer = &top;
bottom_pointer = &bottom;
left_pointer = &left;
right_pointer = &right;

double t[m+2][n+2], tnew[m+1][n+1], diff, difmax;

printf("%d %d %lf\n",m,n, tol);


// start timer: get current time and store it in variable startTIme
gettimeofday(&startTime, NULL);
// initialise temperature array
for (i=0; i <= m+1; i++) {
	for (j=0; j <= n+1; j++) {
		t[i][j] = 30.0;
	}
}

// fix boundary conditions
for (i=1; i <= m; i++) {
	t[i][0] = *left_pointer; // left
	t[i][n+1] = *right_pointer; // right
}
for (j=1; j <= n; j++) {
	t[0][j] = *top_pointer; // top
	t[m+1][j] = *bottom_pointer; // bottom
}

// main loop
iter = 0;
difmax = 1000000.0;
while (difmax > tol) {
	iter++;

	// update temperature for next iteration
	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			tnew[i][j] = (t[i-1][j]+t[i+1][j]+t[i][j-1]+t[i][j+1])/4.0;
		}
	}

	// work out maximum difference between old and new temperatures
	difmax = 0.0;
	for (i=1; i <= m; i++) {
		for (j=1; j <= n; j++) {
			diff = fabs(tnew[i][j]-t[i][j]);
			if (diff > difmax) {
				difmax = diff;
			}
			// copy new to old temperatures
			t[i][j] = tnew[i][j];
		}
	}

}
// stop timer: get current time and store in in variable stopTime
gettimeofday(&stopTime, NULL);


// print results
printf("iter = %d  difmax = %9.11lf", iter, difmax);

/* for (i=0; i <= m+1; i++) {
	printf("\n");
	for (j=0; j <= n+1; j++) {
		printf("%3.5lf ", t[i][j]);
	}
}
printf("\n"); */

// calculate total time by subtracting the startTime from the stopTime (results in microsseconds)
// print the totalTime as a long integer (%ld)
totalTime = (stopTime.tv_sec * 1000000 + stopTime.tv_usec) - (startTime.tv_sec * 1000000 + startTime.tv_usec);
printf("\nThe Total time taken: %ld microseconds\n", totalTime);

}