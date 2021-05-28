/*
============================================================================
Filename    : pi.c
Author      : Jonas Blanc, Mélissa Gehring
SCIPER		: 287508, 264265
============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include "utility.h"

#include <omp.h>

double calculate_pi (int num_threads, int samples);

int main (int argc, const char *argv[]) {

    int num_threads, num_samples;
    double pi;

    if (argc != 3) {
		printf("Invalid input! Usage: ./pi <num_threads> <num_samples> \n");
		return 1;
	} else {
        num_threads = atoi(argv[1]);
        num_samples = atoi(argv[2]);
	}

    set_clock();
    pi = calculate_pi (num_threads, num_samples);

    printf("- Using %d threads: pi = %.15g computed in %.4gs.\n", num_threads, pi, elapsed_time());

    return 0;
}


double calculate_pi (int num_threads, int samples) {
    double pi, x, y = 0;
    
    int countIn = 0;

	int chunk = samples/num_threads;

    omp_set_num_threads(num_threads);
    
    // en //
    #pragma omp parallel private(x, y)
    {
		rand_gen gen = init_rand(); 
		double count = 0;
		
		for(int i = 0; i < chunk; ++i){
			x = next_rand(gen);
			y = next_rand(gen);
			double sumSquare = x*x + y*y;

			if(sumSquare <= 1){
				count++;
			}
		}
		
		free_rand(gen);
		
		#pragma omp critical
		countIn += count;
	}
	
    // en Seq
    pi = (4.0 * countIn) / samples;
    return pi;
}
