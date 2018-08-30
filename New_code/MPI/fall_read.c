#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(void){
	FILE *fp;
	int x;
	long int count = 0;	
	fp = fopen("waterfall_1920_2520.raw","r");
	if(fp == NULL) {
		perror("Error in opening file");
		return(-1);
	} 
	do {
		x = fgetc(fp);
		count++;
		if((count % 1920) == 0){
			printf("\n");
		}
		printf("%d",x);
		if( feof(fp) ) {
			break ;
		}
	} while(1);

	fclose(fp);

	fprintf(stderr, "%ld\n",count);
	return(0);
}