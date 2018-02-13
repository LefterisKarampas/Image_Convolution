
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <cuda.h>
#include <time.h>

#define FRACTION_CEILING(numerator, denominator) ((numerator+denominator-1)/denominator)
#define BLOCK_SIZE  24

__global__ void convolution(unsigned char *in_image, unsigned char *out_image, int height, int width, int *cfilter);


void Usage(char *prog_name) {
   fprintf(stderr, "usage: %s -f <filename> -r <rows> -c <columns> -m <max_loops> -rgb \n", prog_name);
}
 
int main(int argc, char **argv) {
	srand(time(NULL));
	int i, j, n, is_RGB, filter[9] = {1, 2, 1, 2, 4, 2, 1, 2, 1}, *cfilter;
	unsigned int width, height;
	float elapsed_time;
	char *input_path, *output_path;
	unsigned char *Grey_image_array, **RGB_image_array, *tempGrey, **tempRGB, *in_image, *out_image;
	cudaEvent_t start, stop;	

	n = 200;
	i = 1;
	height = 2520;
    width = 1920;
    is_RGB = 0;
	while(i < argc){
      if (!strcmp(argv[i], "-f")){
         input_path = argv[i+1];
      }
      else if (!strcmp(argv[i], "-r"))
      {
         height = atoi(argv[i+1]);
      }
      else if ( !strcmp(argv[i], "-c") )
      {
         width = atoi(argv[i+1]);
      }
      else if ( !strcmp(argv[i], "-m") )
      {
         n = atoi(argv[i+1]);
      }
     else if (!strcmp(argv[i], "-rgb")){
     	is_RGB = 1;
         i--;
     }
      else{
         Usage(argv[0]);
         return 0;
      }
      i+=2;
   }

   	if(!is_RGB){
		Grey_image_array = (unsigned char*) malloc(height * width * sizeof(unsigned char));
		for(i=0;i<height*width;i++){
			Grey_image_array[i] = (unsigned char)rand()%256;
		}
	}
	else{
		RGB_image_array = (unsigned char**) malloc(3 * sizeof(unsigned char*));
		for(i=0;i<3;i++){
			RGB_image_array[i] = (unsigned char*) malloc(height * width * sizeof(unsigned char));
		}
		for (i = 0; i < height * width * 3; i++){
			RGB_image_array[i % 3][i / 3] = (unsigned char)rand()%256;
		}
	}
	
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE,1);
	dim3 dimGrid(FRACTION_CEILING(width, BLOCK_SIZE),FRACTION_CEILING(height, BLOCK_SIZE),1);

	cudaMalloc((void**)&in_image, height * width * sizeof(unsigned char));
	cudaMalloc((void**)&out_image, height * width * sizeof(unsigned char));
	cudaMalloc((void**)&cfilter, 9 * sizeof(int));

	cudaMemcpy(cfilter, filter, 9 * sizeof(int), cudaMemcpyHostToDevice);
	
	//time starts
	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	
	
	if(!is_RGB){
		cudaMalloc((void**)&tempGrey, height * width * sizeof(unsigned char));
		cudaMemcpy(in_image, Grey_image_array, height * width * sizeof(unsigned char), cudaMemcpyHostToDevice);
		for(i = 0; i < n; i++){
			convolution<<<dimGrid,dimBlock>>>(in_image, out_image, height, width, cfilter);
			cudaThreadSynchronize();
			tempGrey = in_image;
			in_image = out_image;
			out_image = tempGrey;
		}
		cudaMemcpy(Grey_image_array, out_image, height * width * sizeof(unsigned char), cudaMemcpyDeviceToHost);
		cudaFree(tempGrey);
	}else{
		tempRGB = (unsigned char**)malloc(3 * sizeof(unsigned char*));
		for(i = 0; i < 3; i++)
			cudaMalloc((void**)&tempRGB[i], height * width * sizeof(unsigned char));
		for(i = 0; i < 3; i++){
			cudaMemcpy(in_image, RGB_image_array[i], height * width * sizeof(unsigned char), cudaMemcpyHostToDevice);
			for(j = 0; j < n; j++){
				convolution<<<dimGrid,dimBlock>>>(in_image, out_image, height, width, cfilter);
				cudaThreadSynchronize();
				tempRGB[i] = in_image;
				in_image = out_image;
				out_image = tempRGB[i];
			}
			cudaMemcpy(RGB_image_array[i], out_image, height * width * sizeof(unsigned char), cudaMemcpyDeviceToHost);
		}
		for(i = 0; i < 3; i++)
			cudaFree(tempRGB[i]);
		free(tempRGB);
	}

	//time finishes
	cudaEventCreate(&stop);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	
	cudaEventElapsedTime(&elapsed_time, start, stop);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	printf("Finished after %f sec\n", elapsed_time/1000);
	

	cudaFree(in_image);
	cudaFree(out_image);
	cudaFree(cfilter);

	//bye
	exit(EXIT_SUCCESS);
}



/*
 * Function for gpu which apply the filter
 */
__global__ void convolution(unsigned char *in_image, unsigned char *out_image, int height, int width, int *cfilter ) {
	
	const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int y = blockDim.y * blockIdx.y + threadIdx.y;

	if (y >= height || x >= width){
    	   return;
	}
    	

	int i, j, s = 1, y_idx, x_idx, sum = 0;

	for (i = -s; i <= s; i++) {
		for ( j = -s; j <= s; j++) {
			y_idx = y + i;
			x_idx = x + j;
			if (y_idx >= height || y_idx < 0 || x_idx >= width || x_idx < 0) {
				y_idx = y;
				x_idx = x;
			}
			sum += in_image[width*(y_idx)+(x_idx)] * cfilter[3*(i+1)+(j+1)];
		}	
	}

	out_image[width*y+x] =(unsigned char)((float)sum/16);	
}
