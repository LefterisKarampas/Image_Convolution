#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <cuda.h>


#define FRACTION_CEILING(numerator, denominator) ((numerator+denominator-1)/(denominator))


__global__ void convolution(int *in,int *out,float *h,int N,int num_elements){
	int index = threadIdx.x + blockDim.x * blockIdx.x;
	if((threadIdx.x == 0) || (threadIdx.x == blockDim.x -1) || (blockIdx.x == 0) || (blockIdx.x == N-1)){
		out[index] = 255;
	}
	else{
		int North = index - blockDim.x;
		int South = index + blockDim.x;
		int East = index +num_elements;
		int West = index -num_elements;
		int NE = North +num_elements;
		int NW = North -num_elements;
		int SE = South +num_elements;
		int SW = South -num_elements;
		out[index] = in[North]*h[1] + in[South]*h[7] + in[East]*h[5] + in[West]*h[3] + in[index]*h[4]+
			in[NE]*h[2]+in[NW]*h[0] + in[SE]*h[8]+in[SW]*h[6];

		if(out[index] > 255){
			out[index] = 255;
		}
		else if(out[index] < 0){
			out[index] = 0;
		}
	}
	__syncthreads();
}


void Usage(char *prog_name) {
   fprintf(stderr, "usage: %s -f <filename> -r <rows> -c <columns> -m <max_loops> -rgb\n", prog_name);
} 

int main(int argc,char **argv){

	int i = 1;
	int N = 2048;
	int M = 1024;
	int max_loops = 200;
	int num_elements = 1;
	char * filename;
	while(i<argc){
		if (!strcmp(argv[i], "-f")){
			filename = argv[i+1];
		}
		else if (!strcmp(argv[i], "-r"))
		{
			N = atoi(argv[i+1]);
		}
		else if ( !strcmp(argv[i], "-c") )
		{
			M = atoi(argv[i+1]);
		}
		else if ( !strcmp(argv[i], "-m") )
		{
			max_loops = atoi(argv[i+1]);
		}
		else if (!strcmp(argv[i], "-rgb")){
			num_elements = 3;
			i--;
		}
		else{
			Usage(argv[0]);
			return 0;
		}
		i+=2;
   	}	

	srand(time(NULL));
	float h[9] = {0,1,3,0,1,0,0,0,0};
	int size = N*M*sizeof(int)*num_elements;
	int *Image = (int *)malloc(size);
	int *Out = (int *)malloc(size);
	int *out;
	int *in;
	float *h_d;
	float elapsed_time;
	cudaEvent_t start, stop;


	dim3 NumberOfThreads(M*num_elements);
	dim3 NumberOfBlocks(N);

	cudaMalloc((void **)&in,size);
	cudaMalloc((void **)&out,size);
	cudaMalloc((void**)&h_d,9*sizeof(float));

	for(int i=0;i<N*M*num_elements;i++){
		Image[i] = rand() % 256;
	}

	cudaMemcpy(h_d,h,9*sizeof(float),cudaMemcpyHostToDevice);

	cudaEventCreate(&start);
	cudaEventRecord(start,0);

	int loop = 0;
	do{
		if(loop % 2 == 0){
			cudaMemcpy(in,Image,size,cudaMemcpyHostToDevice);
			convolution<<<NumberOfBlocks,NumberOfThreads>>>(in,out,h_d,N,num_elements);
			cudaMemcpy(Out,out,size,cudaMemcpyDeviceToHost);
		}
		else{
			cudaMemcpy(in,Out,size,cudaMemcpyHostToDevice);	
			convolution<<<NumberOfBlocks,NumberOfThreads>>>(in,out,h_d,N,num_elements);
			cudaMemcpy(Image,out,size,cudaMemcpyDeviceToHost);
		}
		loop++;
	}while(loop < max_loops && memcmp(Image,Out,size) != 0);

	/*for(int i=0;i<size;i++){
		printf("%d - %d\n",Image[i],Out[i]);
	}*/
	cudaEventCreate(&stop);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	
	cudaEventElapsedTime(&elapsed_time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time %f sec with loop %d\n",elapsed_time/1000,loop);
	free(Image);
	free(Out);
	cudaFree(h_d);
	cudaFree(out);
	cudaFree(in);
	return 0;
}