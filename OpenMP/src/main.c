#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */ 
#include <stdlib.h>
#include <time.h>
#include <math.h>


#ifdef _OPENMP
#include <omp.h>
#endif


static inline int convolution(short int **Table,short int **Final,int i,int j,float h[3][3],int num_elements){
   Final[i][j] = (short int)(h[0][0] * Table[i-1][j-num_elements]) + (short int)(h[0][1] * Table[i-1][j]) + (short int)(h[0][2]*Table[i-1][j+num_elements]) +
      (short int)(Table[i][j-num_elements] * h[1][0]) + (short int)(h[1][1]*Table[i][j]) + (short int)(h[1][2] * Table[i][j+num_elements]) + 
      (short int)(h[2][0]*Table[i+1][j-num_elements])+(short int)(h[2][1] *Table[i+1][j]) + (short int)(h[2][2] *Table[i+1][j+num_elements]);
   return (Final[i][j] == Table[i][j]);
}

void Usage(char *prog_name) {
   fprintf(stderr, "usage: %s -f <filename> -r <rows> -c <columns> -m <max_loops> -rgb -t <threads>\n", prog_name);
}  /* Usage */


int main(int argc,char **argv) {
   srand(time(NULL));
   int comm_sz;
   int my_rank;
   int i;

   //Create Filter
   short int k[3][3] = {{1,0,0},{0,2,0},{0,0,1}};
   float h[3][3];
   int sum;
   int max;
   for(int i=0;i<3;i++){
      sum = 0;
      for(int j=0;j<3;j++){
         sum += k[i][j];
      }
      if(i == 0){
         max = sum;
      }
      else if(sum > max){
         max = sum;
      }
   }
   for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
         h[i][j] = k[i][j]/(float)max;
      }
   }


   /* Start up MPI */
   MPI_Init(NULL, NULL);

   /* Get the number of processes */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   int thr = 7;
   int N = 2520;
   int M = 1920;
   int max_loops = 200;
   char * filename = NULL;
   int num_elements = 1;
   i = 1;
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
      else if (!strcmp(argv[i], "-t"))
      {
         thr = atoi(argv[i+1]);
      }
     else if (!strcmp(argv[i], "-rgb")){
         num_elements = 3;
         i--;
     }
      else{
         if(my_rank == 0){
            Usage(argv[0]);
         }
         MPI_Finalize();
         return 0;
      }
      i+=2;
   }

   if(my_rank == 0){
      fprintf(stderr,"Threads: %d\n",thr);
      fprintf(stderr,"N: %d\n",N);
      fprintf(stderr,"M: %d\n",M);
      fprintf(stderr,"max_loops: %d\n",max_loops);
      if(filename != NULL){
         fprintf(stderr,"filename: %s\n",filename);
      }
   }

   int line_div, col_div;
   double processors_sqrt = sqrt(comm_sz);
   double div;
   int last_div_ok;

   if(processors_sqrt == floor(processors_sqrt)){
      line_div = processors_sqrt;
      col_div = processors_sqrt;
   }
   else{
      if(my_rank == 0)
         fprintf(stderr,"Failed share the data\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
   }


   int rows_per_block = N / line_div;
   int cols_per_block = M / col_div;
   int blocks_per_row = N / rows_per_block;
   int blocks_per_col = M / cols_per_block;

   //Define cols datatype
   MPI_Datatype cols_type;
   MPI_Type_vector(rows_per_block,num_elements,(cols_per_block+2)*num_elements,MPI_SHORT,&cols_type);
   MPI_Type_commit(&cols_type);


   /* Get my rank among all the processes */
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


   if (my_rank == 0) {
      printf("line div : %d\n", line_div);
      printf("col div : %d\n", col_div);
      printf("rows_per_block: %d\n", rows_per_block);
      printf("cols_per_block: %d\n", cols_per_block);
      printf("blocks_per_row: %d\n",blocks_per_row);
      printf("blocks_per_col: %d\n",blocks_per_col);

   }

   short int ** Table = (short int **)malloc(sizeof(short int *)*(rows_per_block+2)); 
   short int ** Final = (short int **)malloc(sizeof(short int *)*(rows_per_block+2));
   for(int i=0;i<rows_per_block+2;i++){
      Table[i] = (short int *)malloc(sizeof(short int)*(cols_per_block+2)*num_elements);
      Final[i] = (short int *)malloc(sizeof(short int)*(cols_per_block+2)*num_elements);
      if((i == 0) || (i == (rows_per_block+1))){
       for(int j=0;j<(cols_per_block+2)*num_elements;j++){
          Table[i][j] = 255;
       }
       continue;
      }
      for(int j=1;j<(cols_per_block+1)*num_elements;j++){
         Table[i][j] = rand() % 256;
      }
   }


   //Neighbors

   MPI_Comm comm;
   int num_dims = 2;
   int reorder = 1;
   int dim[2] = {line_div,line_div};
   int period[2] = {0,0};
   int coord[2];
   int neigh[2];
   MPI_Cart_create(MPI_COMM_WORLD, num_dims, dim, period, reorder, &comm); 
   MPI_Cart_coords(comm, my_rank, num_dims, coord);


   int North,South,East,West,NW,NE,SW,SE;
   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1];
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      North = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      North = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &North);
   }

   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1];
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      South = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      South = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &South);
   }

   neigh[0] = coord[0];
   neigh[1] = coord[1] + 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      East = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      East = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &East);
   }

   neigh[0] = coord[0];
   neigh[1] = coord[1] - 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      West = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      West = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &West);
   }

   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1] + 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      NE = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      NE = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &NE);
   }

   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1] - 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      NW = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      NW = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &NW);
   }

   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1] + 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      SE = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      SE = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &SE);
   }

   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1] - 1;
   if(neigh[0] < 0 || neigh[0] >= blocks_per_col){
      SW = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0 || neigh[1] >= blocks_per_row){
      SW = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &SW);
   }
   /*
   printf("Rank: %d -> (%d,%d)\n \tNorth: %d, South: %d, West: %d, East: %d\n\t\
       North West: %d, North East: %d, South West: %d, South East: %d\n",
      my_rank,coord[0],coord[1],North,South,West,East,NW,NE,SW,SE);
   */

   //Create messages
   MPI_Request send_request[8];
   int tag = 11;
   //rows
   MPI_Send_init(&Table[1][num_elements],cols_per_block*num_elements, MPI_SHORT, North, tag, comm, &send_request[0]);
   MPI_Send_init(&Table[rows_per_block][num_elements],cols_per_block*num_elements, MPI_SHORT,South,tag, comm, &send_request[1]);
   //cols
   MPI_Send_init(&Table[1][num_elements],1,cols_type,West,tag,comm,&send_request[2]);
   MPI_Send_init(&Table[1][cols_per_block*num_elements],1,cols_type,East,tag,comm,&send_request[3]);
   //corners
   MPI_Send_init(&Table[1][num_elements],num_elements, MPI_SHORT,NW, tag, comm, &send_request[4]);
   MPI_Send_init(&Table[1][cols_per_block*num_elements],num_elements,MPI_SHORT,NE,tag,comm,&send_request[5]);
   MPI_Send_init(&Table[rows_per_block][num_elements],num_elements,MPI_SHORT,SW,tag,comm,&send_request[6]);
   MPI_Send_init(&Table[rows_per_block][cols_per_block*num_elements],num_elements,MPI_SHORT,SE,tag, comm, &send_request[7]);

   MPI_Request receive_request[8];
   //rows
   MPI_Recv_init(&Table[0][num_elements],cols_per_block*num_elements, MPI_SHORT, North, tag, comm, &receive_request[0]);
   MPI_Recv_init(&Table[rows_per_block+1][num_elements], cols_per_block*num_elements, MPI_SHORT, South, tag, comm, &receive_request[1]);
   //cols
   MPI_Recv_init(&Table[1][0],1, cols_type, West, tag, comm, &receive_request[2]);
   MPI_Recv_init(&Table[1][(cols_per_block+1)*num_elements],1, cols_type, East, tag, comm, &receive_request[3]);
   //corners
   MPI_Recv_init(&Table[0][0],num_elements,MPI_SHORT,NW,tag,comm,&receive_request[4]);
   MPI_Recv_init(&Table[0][(cols_per_block+1)*num_elements],num_elements,MPI_SHORT,NE,tag,comm,&receive_request[5]);
   MPI_Recv_init(&Table[rows_per_block+1][0],num_elements,MPI_SHORT,SW,tag,comm,&receive_request[6]);
   MPI_Recv_init(&Table[rows_per_block+1][(cols_per_block+1)*num_elements],num_elements,MPI_SHORT,SE,tag,comm,&receive_request[7]);


   MPI_Status status[8];
   int loop = 0;
   int changes = 0;
   //int lchanges = 0;
   int sum_changes,j;
   double start, finish;
   int lchanges=0;
   omp_set_num_threads(thr);
   start = MPI_Wtime();
   while(loop < max_loops){
      loop++;
      //8x ISend
      MPI_Startall(8, send_request);
      //8x IRecv
      MPI_Startall(8, receive_request);
      //Do for our table
      if (changes>0)
        lchanges = 1;
  #pragma omp parallel for reduction(+:changes) private(lchanges)
  //{
      for(int i=2;i<rows_per_block;i++){
     //#pragma omp for schedule(dynamic) reduction(+:changes) private(lchanges)
  	for(int j=2*num_elements;j<cols_per_block*num_elements;j++){
          if(convolution(Table, Final,i,j,h,num_elements) && !lchanges){
               lchanges++;
               changes+=lchanges;
            }
         }
      }
  //}
      MPI_Waitall(8, receive_request,status);
      //do the job for receive
      //First row
   if (changes>0)
        lchanges = 1;
  #pragma omp parallel for schedule(dynamic) reduction(+:changes) private(lchanges)
    for(j=num_elements;j<(cols_per_block+1)*num_elements;j++){
	//printf("thread %d and j %d\n",omp_get_thread_num(),j);
         if(convolution(Table, Final,1,j,h,num_elements) && !lchanges){
	         lchanges++;
            changes+=lchanges;
         }
    }
    if (changes>0)
	lchanges = 1;
  #pragma omp parallel for schedule(dynamic) reduction(+:changes) private(lchanges)
      //Last row
      for(int j=num_elements;j<(cols_per_block+1)*num_elements;j++){
         if(convolution(Table, Final,rows_per_block,j,h,num_elements) && !lchanges){
            lchanges++;
            changes+=lchanges;
         }
      }
    if (changes>0)
        lchanges = 1;
  #pragma omp parallel for schedule(dynamic) reduction(+:changes) private(lchanges)
      //First col and last col for each middle row
      for(int j=2;j<rows_per_block;j++){
         for(int k=0;k<num_elements;k++){
            if(convolution(Table,Final,j,num_elements+k,h,num_elements) && !lchanges){
               lchanges++;
               changes+=lchanges;
            }

           if (convolution(Table,Final,j,cols_per_block*num_elements+k,h,num_elements) && !lchanges){
                   lchanges++;
                   changes+=lchanges;
           }
         }
     }
      short int ** temp;
      temp = Table;
      Table = Final;
      Final = temp;
      MPI_Waitall(8,send_request,status);
      //Reduce all changes
      if(loop % 10 == 0){
         MPI_Allreduce(&changes,&sum_changes, 1, MPI_SHORT, MPI_SUM, comm);
         changes = lchanges = 0;
         if(sum_changes == 0){
            if(my_rank == 0)
              printf("No changes in loop %d!\n",loop);
           // break;
            //loop = max_loops;
         }
         /*else if(loop > 100){
            h[0][1] = 0;
         }*/
      }
   }
   finish = MPI_Wtime();

   double local_elapsed, elapsed;
   local_elapsed = finish - start;
   MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   if(my_rank == 0){
      printf("Time elapsed: %f seconds\n", elapsed);
   }
   /* Shut down MPI */
   for(int i=0;i<line_div;i++){
      free(Table[i]);
      free(Final[i]);
   }
   free(Table);
   free(Final);
   MPI_Finalize();

   return 0;
}  /* main */
