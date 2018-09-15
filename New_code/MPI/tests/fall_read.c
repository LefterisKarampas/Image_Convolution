#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc,char *argv[]){
	FILE *fp;
	int x;
	long int count = -1;	
	fp = fopen(argv[2],"r");
	if(fp == NULL) {
		perror("Error in opening file");
		return(-1);
	} 
	do {
		x = fgetc(fp);
		count++;
		if((count % atoi(argv[1])) == 0){
			printf("\n");
		}
		printf("%c",x);
		if( feof(fp) ) {
			break ;
		}
	} while(1);

	fclose(fp);
	return(0);
}