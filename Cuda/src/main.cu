#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>



__global__ void convolution(int *in,int *out,float *h,int N){
	int index = threadIdx.x + blockDim.x * blockIdx.x;
	if((threadIdx.x == 0) || (threadIdx.x == blockDim.x -1) || (blockIdx.x == 0) || (blockIdx.x == N-1)){
		out[index] = 255;
	}
	else{
		int North = index - blockDim.x;
		int South = index + blockDim.x;
		int East = index +1;
		int West = index -1;
		int NE = North +1;
		int NW = North -1;
		int SE = South +1;
		int SW = South -1;
		out[index] = in[North]*h[1] + in[South]*h[7] + in[East]*h[5] + in[West]*h[3] + in[index]*h[4]+
			in[NE]*h[2]+in[NW]*h[0] + in[SE]*h[8]+in[SW]*h[6];

		if(out[index] > 255){
			out[index] = 255;
		}
		else if(out[index] < 0){
			out[index] = 0;
		}

	}
}

#define N 2048
#define M 1024


int main(int argc,char **argv){
	srand(time(NULL));
	float h[9] = {0,0,0,0,1,0,0,0,0};
	int size = N*M*sizeof(int);
	int *Image = (int *)malloc(size);
	int *Out = (int *)malloc(size);
	int *out;
	int *in;
	float *h_d;
	cudaMalloc((void **)&in,size);
	cudaMalloc((void **)&out,size);
	cudaMalloc((void**)&h_d,9*sizeof(float));
	for(int i=0;i<N*M;i++){
		Image[i] = rand() % 256;
	}
	cudaMemcpy(h_d,h,9*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(in,Image,size,cudaMemcpyHostToDevice);


	convolution<<<N,M>>>(in,out,h_d,N);
	cudaMemcpy(Out,out,size,cudaMemcpyDeviceToHost);

	/*int loop = 1;
	do{
		cudaMemcpy(Out,out,size,cudaMemcpyDeviceToHost);
		cudaMemcpy(in,Out,size,cudaMemcpyHostToDevice);
		convolution<<<N,M>>>(in,out,h_d,N);
		loop++;
	}while(loop < 10);
	*/

	for(int i=0;i<N;i++){
		for(int j=0;j<M;j++){
			printf("%3d ",Out[i*M+j]);
		}
		printf("\n");
	}

	free(Image);
	free(Out);
	cudaFree(h_d);
	cudaFree(out);
	cudaFree(in);
	return 0;
}