#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */ 
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 2520
#define M 1920

static inline void convolution(short int **Table,short int **Final,int i,int j,float h[3][3]){
   Final[i][j] = (short int)(h[0][0] * Table[i-1][j-1]) + (short int)(h[0][1] * Table[i-1][j]) + (short int)(h[0][2]*Table[i-1][j+1]) +
      (short int)(Table[i][j-1] * h[1][0]) + (short int)(h[1][1]*Table[i][j]) + (short int)(h[1][2] * Table[i][j+1]) + 
      (short int)(h[2][0]*Table[i+1][j-1])+(short int)(h[2][1] *Table[i+1][j]) + (short int)(h[2][2] *Table[i+1][j+1]);
}

int main(void) {
   srand(time(NULL));
   int comm_sz;          
   int my_rank;
   int i;

   //Create Filter
   short int k[3][3] = {{0,1,0},{0,1,0},{0,0,0}};
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
   MPI_Type_vector(rows_per_block,1,0,MPI_SHORT,&cols_type);
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
      Table[i] = (short int *)malloc(sizeof(short int)*(cols_per_block+2));
      Final[i] = (short int *)malloc(sizeof(short int)*(cols_per_block+2));
      if((i == 0) || (i == (rows_per_block+1))){
         /*for(int j=0;j<cols_per_block+2;j++){
            Table[i][j] = 255;
         }*/
         continue;
      }
      for(int j=1;j<cols_per_block+1;j++){
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
   MPI_Send_init(&Table[1][1],cols_per_block, MPI_SHORT, North, tag, comm, &send_request[0]);
   MPI_Send_init(&Table[rows_per_block][1],cols_per_block, MPI_SHORT,South,tag, comm, &send_request[1]);
   //cols
   MPI_Send_init(&Table[1][1],1,cols_type,West,tag,comm,&send_request[2]);
   MPI_Send_init(&Table[1][cols_per_block],1,cols_type,East,tag,comm,&send_request[3]);
   //corners
   MPI_Send_init(&Table[1][1],1, MPI_SHORT,NW, tag, comm, &send_request[4]);
   MPI_Send_init(&Table[1][cols_per_block],1,MPI_SHORT,NE,tag,comm,&send_request[5]);
   MPI_Send_init(&Table[rows_per_block][1],1,MPI_SHORT,SW,tag,comm,&send_request[6]);
   MPI_Send_init(&Table[rows_per_block][cols_per_block],1,MPI_SHORT,SE,tag, comm, &send_request[7]);

   MPI_Request receive_request[8];
   //rows
   MPI_Recv_init(&Table[0][1],cols_per_block, MPI_SHORT, North, tag, comm, &receive_request[0]);
   MPI_Recv_init(&Table[rows_per_block+1][1], cols_per_block, MPI_SHORT, South, tag, comm, &receive_request[1]);
   //cols
   MPI_Recv_init(&Table[1][0],1, cols_type, West, tag, comm, &receive_request[2]);
   MPI_Recv_init(&Table[1][cols_per_block+1], 1, cols_type, East, tag, comm, &receive_request[3]);
   //corners
   MPI_Recv_init(&Table[0][0],1,MPI_SHORT,NW,tag,comm,&receive_request[4]);
   MPI_Recv_init(&Table[0][cols_per_block+1],1,MPI_SHORT,NE,tag,comm,&receive_request[5]);
   MPI_Recv_init(&Table[rows_per_block+1][0],1,MPI_SHORT,SW,tag,comm,&receive_request[6]);
   MPI_Recv_init(&Table[rows_per_block+1][cols_per_block+1],1,MPI_SHORT,SE,tag,comm,&receive_request[7]);
  

   MPI_Status status[8];
   int loop = 0;
   int changes = 0;
   int sum_changes;
   double start, finish;

   start = MPI_Wtime();
   while(loop < 100){
      changes = 0;
      loop++;

      //8x ISend
      MPI_Startall(8, send_request);
      //8x IRecv
      MPI_Startall(8, receive_request);
      
      //Do for our table
      for(int i=2;i<rows_per_block;i++){
         for(int j=2;j<cols_per_block;j++){
            convolution(Table, Final,i,j,h);
            if(!changes && Final[i][j] != Table[i][j]){
               changes++;
            }
         }
      }
      MPI_Waitall(8, receive_request,status);

      //do the job for receive
      for(int i=1;i<rows_per_block+1;i++){
         if((i == 1) || (i == rows_per_block)){
            for(int j=1;j<cols_per_block+1;j++){
               convolution(Table, Final,i,j,h);
               if(!changes && Final[i][j] != Table[i][j]){
                  changes++;
               }
            }
         }
         else{
            int j = 1;
            convolution(Table,Final,i,j,h);
            if(!changes && Final[i][j] != Table[i][j]){
               changes++;
            }
            j = cols_per_block;
            convolution(Table,Final,i,j,h);
            if(!changes && Final[i][j] != Table[i][j]){
               changes++;
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
         if(sum_changes == 0){
            if(my_rank == 0)
               printf("No changes in loop %d!\n",loop);
            break;
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
