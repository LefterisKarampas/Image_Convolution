#ifndef _HALLO_POINTS_H_
#define _HALLO_POINTS_H_


typedef struct halo{
	unsigned char * North_East;
	unsigned char * North_West;
	unsigned char * South_East;
	unsigned char * South_West;
	unsigned char * North;
	unsigned char * South;
	unsigned char * West;
	unsigned char * East;
}Halo_points;


void Initialize_Halo(Halo_points *,int,int,int);

void Delete_Halo(Halo_points *);

#endif