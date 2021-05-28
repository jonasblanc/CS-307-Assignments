/*
============================================================================
Filename    : algorithm.c
Author      : Jonas Blanc, Mélissa Gehring
SCIPER		: 287508, 264265
============================================================================
*/
#include <math.h>
#include <stdlib.h>


#define INPUT(I,J) input[(I)*length+(J)]
#define OUTPUT(I,J) output[(I)*length+(J)]

void simulate_parallel_standard(double *input, double *output, int threads, int length, int iterations);
void simulate_parallel_inverted(double *input, double *output, int threads, int length, int iterations);
void simulate_axis_inverted(double *input, double *output, int threads, int length, int iterations);

void simulate_matrix2_hori(double *input, double *output, int threads, int length, int iterations);
void simulate_matrix(double *input, double *output, int threads, int length, int iterations);

void simulate_unrolling_hori(double *input, double *output, int threads, int length, int iterations);
void simulate_unrolling_verti(double *input, double *output, int threads, int length, int iterations);


//====================================================== What we tried before knowing it was not allowed ======================================================//

/*void simulate_divide_array_in_4(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_array_in_8(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_array_in_2_hori(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_array_in_2_verti(double *input, double *output, int threads, int length, int iterations);

void simulate_grow_whit_it(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_4_through_it(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_4_through_and_grow_it(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_8_through_it(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_8_through_it_dynamic(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_8_through_and_grow_it(double *input, double *output, int threads, int length, int iterations);
void simulate_divide_8_through_and_grow_it_dynamic(double *input, double *output, int threads, int length, int iterations);
*/

void simulate(double *input, double *output, int threads, int length, int iterations)
{
	/*
	 * Playing with parallelization and loop order
	 */
    //simulate_parallel_standard(input, output, threads, length,iterations);
    //simulate_parallel_inverted(input, output, threads, length,iterations);
    //simulate_axis_inverted(input, output, threads, length,iterations);
    
    /*
     * Attempt to optimize locality
     */
    //simulate_matrix2_hori(input, output, threads, length,iterations);
    //simulate_matrix(input, output, threads, length,iterations);
    
    /*
    * Unrolling principle: do more computation per loop to reduce number of loops
    */
    simulate_unrolling_hori(input, output, threads, length,iterations);
    //simulate_unrolling_verti(input, output, threads, length,iterations);


    //====================================================== Not allowed ======================================================//
    /*
    * Symetry principle: Calculate only once for every symetry axes
    */
    //simulate_divide_array_in_4(input, output, threads, length,iterations);
    //simulate_divide_array_in_8(input, output, threads, length,iterations);
    //simulate_divide_array_in_2_hori(input, output, threads, length,iterations);
    //simulate_divide_array_in_2_verti(input, output, threads, length,iterations);

    /*
     * Grow with iteration principle: We  don't need to iterate over all the array since heat only propagate to next cell every iteration
     */
    //simulate_grow_whit_it(input, output, threads, length,iterations);

    /*
     * No extra write principle: using symetry principle, no need to write to all the array at each iteration, copy to all the array at the end
     */
    //simulate_divide_4_through_it(input, output, threads, length,iterations);
    //simulate_divide_8_through_it(input, output, threads, length,iterations);
    //simulate_divide_8_through_it_dynamic(input, output, threads, length,iterations);

    /*
     * Best algo so far: combine Grow with iteration principle + No extra write principle
     */
    //simulate_divide_4_through_and_grow_it(input, output, threads, length,iterations);
    //simulate_divide_8_through_and_grow_it(input, output, threads, length,iterations);
    //simulate_divide_8_through_and_grow_it_dynamic(input, output, threads, length,iterations);

    //====================================================== Not allowed ======================================================//
}

void simulate_parallel_standard(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		#pragma omp parallel for 
        for(int i=1; i<length-1; i++)
        {
            for(int j=1; j<length-1; j++)
            {
                    if ( ((i == length/2-1) || (i == length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_parallel_inverted(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		for(int i=1; i<length-1; i++)
        {
			#pragma omp parallel for
            for(int j=1; j<length-1; j++)
            {
				printf("%d\n", omp_get_thread_num());
                    if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_axis_inverted(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		#pragma omp parallel for 
        for(int j=1; j<length-1; j++)
        {
            for(int i=1; i<length-1; i++)
            {
                    if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_matrix2_hori(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		#pragma omp parallel
		{ 
			int id = omp_get_thread_num();
			int matrix_size = length / threads;
			int min = id * matrix_size;
			int max = min + matrix_size;
			for(int n = 0; n < threads; n++) {
				for(int i = min; i < max; i++)
				{
					for(int j = n*matrix_size; j < (n+1)*matrix_size; j++)
					{
						if(i != 0 && i != length-1 && j != 0 && j != length -1){
							if ( ((i == length/2-1) || (i == length/2))
								&& ((j == length/2-1) || (j == length/2)) )
								continue;

							OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
										   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
										   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
						}
					}
				}
			}
		}

        temp = input;
        input = output;
        output = temp;
    }
}

int shouldBeCalculate(int i, int j, int lengthMinusOne, int midLengthMinusOne, int midLength){
   
    if(i < 1 || i >= lengthMinusOne || j < 1 || lengthMinusOne <= j){
        return 0;
    }
    if((i == midLengthMinusOne || i== midLength)
            && (j == midLengthMinusOne || j == midLength)){
        return 0;
    }
    return 1;
}

void simulate_matrix(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;

    int lengthMinusOne = length -1;
    int midLength = length / 2;
    int midLengthMinusOne = midLength -1;


    int x_length_block = 16;
    int y_length_block = 25;

    // Number of blocks needed to cover the whole array on the x-axis
    int block_per_row = length / x_length_block;
    if(length % x_length_block != 0){
        block_per_row++;
    }

    // Number of blocks needed to cover the whole array on the y-axis
    int block_per_column = length / y_length_block;
    if(length % y_length_block != 0){
        block_per_column++;
    }

    omp_set_num_threads(threads);


    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        // Iterate along x-axis blocks
        #pragma omp parallel for
        for(int x_block_idx = 0; x_block_idx < block_per_row; x_block_idx++)
        {
            int offset_in_row = x_block_idx * x_length_block; // Offset of current block (ie where it is on x-axis )

            int offset_in_column = 0; // Offset of current block (ie where it is on y-axis )
            // Iterate along y-axis blocks
            for(int y_block_idx = 0; y_block_idx < block_per_column; y_block_idx++)
            {
                int j = offset_in_row; // Position of current element in matrix on x-axis
                // Iterate over elements in current block over x-axis
                for(int x_pos_in_cur_block= 0; x_pos_in_cur_block < x_length_block; ++x_pos_in_cur_block)  
                {
                    int i = offset_in_column; // Position of current element in matrix on y-axis
                    // Iterate over elements in current block over y-axis
                    for(int y_pos_in_cur_block = 0; y_pos_in_cur_block < y_length_block; ++y_pos_in_cur_block)
                    {   
                        // If the current element is in the matrix and not in the middle (source) compute it's output
                        if(shouldBeCalculate(i,j, lengthMinusOne, midLengthMinusOne, midLength) == 1){
                            OUTPUT(i,j) = 
                                (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                        }
                        i++;
                    }
                    j++;
                }
                offset_in_column += y_length_block;
            }
            offset_in_row += x_length_block;     
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_unrolling_verti(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		#pragma omp parallel for 
        for(int i=1; i<length-1; i+=2)
        {
            for(int j=1; j<length-1; j++)
            {
                    if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    ++i;
                    if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;
                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    i--;
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_unrolling_hori(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);
    
    for(int n=0; n < iterations; n++)
    {
		#pragma omp parallel for 
        for(int i=1; i<length-1; i++)
        {
            for(int j=1; j<length-1; j+=2)
            {
                   if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    ++j;
                    if ( ((i == length/2-1) || (i== length/2))
                        && ((j == length/2-1) || (j == length/2)) )
                        continue;
                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;

                    j--;
                   
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

//====================================================== What we tried before knowing it was not allowed ======================================================//
/*
void simulate_divide_array_in_4(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<midLength; i++)
        {
            for(int j=1; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;
                    OUTPUT(lengthMinusOne-i,j) = result;
                    OUTPUT(i,lengthMinusOne-j) = result;
                    OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
            }
        }
        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_divide_array_in_8(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<midLength; i++)
        {
            for(int j=i; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    if(i == j){
                        OUTPUT(i,j) = result;
                        OUTPUT(lengthMinusOne-i,j) = result;
                        OUTPUT(i,lengthMinusOne-j) = result;
                        OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                    }else{
                        OUTPUT(i,j) = result;
                        OUTPUT(j,i) = result;
                        OUTPUT(lengthMinusOne-i,j) = result;
                        OUTPUT(j, lengthMinusOne-i) = result;
                        OUTPUT(i,lengthMinusOne-j) = result;
                        OUTPUT(lengthMinusOne-j,i) = result;
                        OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                        OUTPUT(lengthMinusOne-j,lengthMinusOne-i) = result;
                    }
            }
        }
        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_divide_array_in_2_hori(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<lengthMinusOne; i++)
        {
            for(int j=1; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;
                    OUTPUT(i,lengthMinusOne - j) = result;
            }
        }
        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_divide_array_in_2_verti(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<midLength; i++)
        {
            for(int j=1; j<lengthMinusOne; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;
                    OUTPUT(lengthMinusOne-i,j) = result;
            }
        }
        temp = input;
        input = output;
        output = temp;
    }
}


void simulate_grow_whit_it(double *input, double *output, int threads, int length, int iterations)
{
    double *temp;
    
    omp_set_num_threads(threads);

    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    
    int bound_start = midLengthMinusOne;
    int bound_end = midLength + 1;

    for(int n=0; n < iterations; n++)
    {
        if(bound_start > 1){
             bound_start--;
        }
        if(bound_end < lengthMinusOne){
             bound_end++;
        }
		#pragma omp parallel for 
        for(int i=bound_start; i<bound_end; i++)
        {
            for(int j=bound_start; j<bound_end; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;

                    OUTPUT(i,j) = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
            }
        }

        temp = input;
        input = output;
        output = temp;
    }
}

void simulate_divide_4_through_it(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<midLength; i++)
        {
            for(int j=1; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    
                    OUTPUT(i,j) = result;
                    if(i == midLengthMinusOne){
                         OUTPUT(midLength,j) = result;
                    }else{
                         if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                        }
                    }
            }
        }
        temp = input;
        input = output;
        output = temp;
    }

    #pragma omp parallel for
    for(int i=1; i<midLength; i++)
    {
        for(int j=1; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
            OUTPUT(lengthMinusOne-i,j) = result;
            OUTPUT(i,lengthMinusOne-j) = result;
            OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
        }
    }
}

void simulate_divide_4_through_and_grow_it(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);

    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    
    int bound_start = midLengthMinusOne;

    for(int n=0; n < iterations; n++)
    {
        if(bound_start > 1){
             bound_start--;
        }

        #pragma omp parallel for
        for(int i=bound_start; i<midLength; i++)
        {
            for(int j=bound_start; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    
                    OUTPUT(i,j) = result;
                    if(i == midLengthMinusOne){
                         OUTPUT(midLength,j) = result;
                    }else{
                         if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                        }
                    }
            }
        }
        temp = input;
        input = output;
        output = temp;
    }

    #pragma omp parallel for
    for(int i=bound_start; i<midLength; i++)
    {
        for(int j=bound_start; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
            OUTPUT(lengthMinusOne-i,j) = result;
            OUTPUT(i,lengthMinusOne-j) = result;
            OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
        }
    }
}

void simulate_divide_8_through_it(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for
        for(int i=1; i<midLength; i++)
        {
            for(int j=i; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;

                    int diff = j-i;
                    if(i != j && diff < 3){
                        OUTPUT(j,i) = result;
                    }

                    if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                    }

            }
        }
        temp = input;
        input = output;
        output = temp;
    }


    #pragma omp parallel for
    for(int i=1; i<midLength; i++)
    {
        for(int j=i; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
             if(i == j){
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
            }else{
                OUTPUT(j,i) = result;
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(j, lengthMinusOne-i) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,i) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,lengthMinusOne-i) = result;
            }
            
        }
    }

}

void simulate_divide_8_through_it_dynamic(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;
    // Parallelize this!!
    for(int n=0; n < iterations; n++)
    {
        #pragma omp parallel for schedule(dynamic, 1)
        for(int i=1; i<midLength; i++)
        {
            for(int j=i; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;

                    int diff = j-i;
                    if(i != j && diff < 3){
                        OUTPUT(j,i) = result;
                    }

                    if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                    }

            }
        }
        temp = input;
        input = output;
        output = temp;
    }


    #pragma omp parallel for
    for(int i=1; i<midLength; i++)
    {
        for(int j=i; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
             if(i == j){
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
            }else{
                OUTPUT(j,i) = result;
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(j, lengthMinusOne-i) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,i) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,lengthMinusOne-i) = result;
            }
            
        }
    }

}

void simulate_divide_8_through_and_grow_it(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;

    int bound_start = midLengthMinusOne;

    for(int n=0; n < iterations; n++)
    {
        if(bound_start > 1){
             bound_start--;
        }
        #pragma omp parallel for
        for(int i=bound_start; i<midLength; i++)
        {
            for(int j=i; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;

                    int diff = j-i;
                    if(i != j && diff < 3){
                        OUTPUT(j,i) = result;
                    }

                    if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                    }

            }
        }
        temp = input;
        input = output;
        output = temp;
    }

    #pragma omp parallel for
    for(int i = bound_start; i < midLength; i++)
    {
        for(int j=i; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
             if(i == j){
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
            }else{
                OUTPUT(j,i) = result;
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(j, lengthMinusOne-i) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,i) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,lengthMinusOne-i) = result;
            }
            
        }
    }

}

void simulate_divide_8_through_and_grow_it_dynamic(double *input, double *output, int threads, int length, int iterations)
{
     double *temp;
    
    omp_set_num_threads(threads);
    
    int midLength = length / 2;
    int midLengthMinusOne = midLength-1;
    int lengthMinusOne = length -1;

    int bound_start = midLengthMinusOne;

    for(int n=0; n < iterations; n++)
    {
        if(bound_start > 1){
             bound_start--;
        }
        #pragma omp parallel for schedule(dynamic, 1)
        for(int i=bound_start; i<midLength; i++)
        {
            for(int j=i; j<midLength; j++)
            {
                    if ( ((i == midLengthMinusOne) || (i== midLength))
                        && ((j == midLengthMinusOne) || (j == midLength)) )
                        continue;
                    double result = (INPUT(i-1,j-1) + INPUT(i-1,j) + INPUT(i-1,j+1) +
                                   INPUT(i,j-1)   + INPUT(i,j)   + INPUT(i,j+1)   +
                                   INPUT(i+1,j-1) + INPUT(i+1,j) + INPUT(i+1,j+1) )/9;
                    OUTPUT(i,j) = result;

                    int diff = j-i;
                    if(i != j && diff < 3){
                        OUTPUT(j,i) = result;
                    }

                    if(j == midLengthMinusOne){
                            OUTPUT(i,midLength) = result;
                    }

            }
        }
        temp = input;
        input = output;
        output = temp;
    }

    #pragma omp parallel for
    for(int i = bound_start; i < midLength; i++)
    {
        for(int j=i; j<midLength; j++)
        {
            if ( ((i == midLengthMinusOne) || (i== midLength))
                && ((j == midLengthMinusOne) || (j == midLength)) )
                    continue;
            double result = OUTPUT(i,j);
             if(i == j){
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
            }else{
                OUTPUT(j,i) = result;
                OUTPUT(lengthMinusOne-i,j) = result;
                OUTPUT(j, lengthMinusOne-i) = result;
                OUTPUT(i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,i) = result;
                OUTPUT(lengthMinusOne-i,lengthMinusOne-j) = result;
                OUTPUT(lengthMinusOne-j,lengthMinusOne-i) = result;
            }
            
        }
    }

}
*/
//====================================================== Not allowed ======================================================//
