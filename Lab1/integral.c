/*
============================================================================
Filename    : integral.c
Author      : Jonas Blanc, Mélissa Gehring
SCIPER		: 287508, 264265
============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include "utility.h"
#include "function.c"

#include <omp.h>

double integrate (int num_threads, int samples, int a, int b, double (*f)(double));

int main (int argc, const char *argv[]) {

    int num_threads, num_samples, a, b;
    double integral;

    if (argc != 5) {
		printf("Invalid input! Usage: ./integral <num_threads> <num_samples> <a> <b>\n");
		return 1;
	} else {
        num_threads = atoi(argv[1]);
        num_samples = atoi(argv[2]);
        a = atoi(argv[3]);
        b = atoi(argv[4]);
	}

    set_clock();

    /* You can use your self-defined funtions by replacing identity_f. */
    integral = integrate (num_threads, num_samples, a, b, identity_f);

    printf("- Using %d threads: integral on [%d,%d] = %.15g computed in %.4gs.\n", num_threads, a, b, integral, elapsed_time());

    return 0;
}


double integrate (int num_threads, int samples, int a, int b, double (*f)(double)) {
    double integral, x = 0;
    double surface = 0;
    
    int interval = b-a;
    
    int chunk = samples/num_threads;

    omp_set_num_threads(num_threads);

    // en //
    #pragma omp parallel private(x)
    {		
		rand_gen gen = init_rand();
		double area = 0; 
		
		for(int i = 0; i < chunk; ++i){
			x = a + next_rand(gen)*interval; //a <= x <= b
			//area += f(x) * interval; 
			area += f(x);
		}
		area *= interval;
		
		free_rand(gen);
		
		#pragma omp critical
		surface += area;
	}
	
    // en Seq
    integral = surface/samples;
    return integral;
}
