#include "../include/halo_points.h"
#include <stdlib.h>



void Initialize_Halo(Halo_points *Halo,int rows,int columns,int num_elements){
	Halo->North = (unsigned char *)malloc(sizeof(unsigned char)*columns*num_elements);
	Halo->South = (unsigned char *)malloc(sizeof(unsigned char)*columns*num_elements);
	Halo->East = (unsigned char *)malloc(sizeof(unsigned char)*rows*num_elements);
	Halo->West = (unsigned char *)malloc(sizeof(unsigned char)*rows*num_elements);

	Halo->North_East = (unsigned char *)malloc(sizeof(unsigned char)*num_elements);
	Halo->North_West = (unsigned char *)malloc(sizeof(unsigned char)*num_elements);
	Halo->South_East = (unsigned char *)malloc(sizeof(unsigned char)*num_elements);
	Halo->South_West = (unsigned char *)malloc(sizeof(unsigned char)*num_elements);
}

void Delete_Halo(Halo_points *Halo){
	free(Halo->North);
	free(Halo->South);
	free(Halo->East);
	free(Halo->West);
	free(Halo->North_East);
	free(Halo->North_West);
	free(Halo->South_East);
	free(Halo->South_West);
}

