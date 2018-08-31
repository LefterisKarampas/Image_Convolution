#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */ 
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../include/halo_points.h"

static inline int convolution(unsigned char **Table,unsigned char **Final,int i,int j,float h[3][3],int num_elements){
   Final[i][j] = (unsigned char) ceil(h[0][0] * Table[i-1][j-num_elements] + h[0][1] * Table[i-1][j] + h[0][2]*Table[i-1][j+num_elements] +
     Table[i][j-num_elements] * h[1][0] + (h[1][1]*Table[i][j]) + (h[1][2] * Table[i][j+num_elements]) + 
      (h[2][0]*Table[i+1][j-num_elements])+(h[2][1] *Table[i+1][j]) + (h[2][2] *Table[i+1][j+num_elements])) % 256;
   return (Final[i][j] == Table[i][j]);
}

void Usage(char *prog_name) {
   fprintf(stderr, "usage: %s -f <filename> -r <rows> -c <columns> -m <max_loops> -o <output_file>\n", prog_name);
}  /* Usage */


int main(int argc,char **argv) {
   srand(time(NULL));
   int comm_sz;          
   int my_rank;
   int i;

   //Create Filter
   unsigned char k[3][3] = {{0,0,0},{0,1,0},{0,0,0}};
   float h[3][3];
   int sum = 0;
   for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
         sum += (int)k[i][j];
      }
   }
   for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
         h[i][j] = k[i][j]/(float)sum;
      }
   }
   

   /* Start up MPI */
   MPI_Init(NULL, NULL); 


   /* Get the number of processes */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   /* Read Program Arguments*/
   int N = 2520;
   int M = 1920;
   int max_loops = 200;
   char * filename = NULL;
   char * output = NULL;
   int num_elements = 1;
   i=1;
   while(i<argc){
      if (!strcmp(argv[i], "-f")){
         filename = argv[i+1];
      }
      else if (!strcmp(argv[i], "-o")){
         output = argv[i+1];
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
      else if ( !strcmp(argv[i], "-rgb") )
      {
         num_elements = 3;
         if(my_rank == 0){
            printf("num_elements %d\n",num_elements);
         }
         i--;
      }
      else{
         if(my_rank == 0){
            Usage(argv[0]);
         }
         MPI_Finalize();
         return 1;
      }
      i+=2;
   }

   if(filename == NULL){
      if(my_rank == 0){
         Usage(argv[0]);
      }
      MPI_Finalize();
      return 2;
   }

   if(output == NULL){
      if(my_rank == 0){
         Usage(argv[0]);
      }
      MPI_Finalize();
      return 3;
   }

   if(my_rank == 0){
      fprintf(stderr,"N: %d\n",N);
      fprintf(stderr,"M: %d\n",M);
      fprintf(stderr,"max_loops: %d\n",max_loops);
      if(filename != NULL){
         fprintf(stderr,"filename: %s\n",filename);
      }
   }



   int process_row, process_col;
   double processors_sqrt = sqrt(comm_sz);

   if(processors_sqrt == floor(processors_sqrt)){
      process_row = processors_sqrt;
      process_col = processors_sqrt;
   }
   else{
      if(my_rank == 0)
         fprintf(stderr,"Failed share the data\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
   }


   int rows_per_block = ceil(N/(float)process_row);
   int cols_per_block = ceil(M /(float)process_col);




   if (my_rank == 0) {
      fprintf(stderr,"rows_per_block: %d\n",rows_per_block);
      fprintf(stderr,"cols_per_block: %d\n",cols_per_block*num_elements);
      fprintf(stderr,"process_row: %d\n",process_row);
      fprintf(stderr,"process_col: %d\n",process_col);
   }
   

   //Allocate Process Tables and Halo Points
   unsigned char ** Table = (unsigned char **)malloc(sizeof(unsigned char *)*rows_per_block);
   Table[0] = (unsigned char *)malloc(sizeof(unsigned char)*rows_per_block*cols_per_block*num_elements); 
   unsigned char ** Final = (unsigned char **)malloc(sizeof(unsigned char *)*rows_per_block);
   Final[0] = (unsigned char *)malloc(sizeof(unsigned char)*rows_per_block*cols_per_block*num_elements);
  
  for(int i=0;i<rows_per_block;i++){
      Table[i]=&(Table[0][i*cols_per_block*num_elements]);
      Final[i]=&(Final[0][i*cols_per_block*num_elements]);
   }


   Halo_points * halo_p = (Halo_points *)malloc(sizeof(Halo_points));
   Initialize_Halo(halo_p,rows_per_block,cols_per_block,num_elements);
   

   //Neighbors
   MPI_Comm comm;
   int num_dims = 2;
   int reorder = 1;
   int dim[2] = {process_row,process_col}; 
   int period[2] = {0,0};
   int coord[2];
   int neigh[2];
   MPI_Cart_create(MPI_COMM_WORLD, num_dims, dim, period, reorder, &comm); 
   MPI_Cart_coords(comm, my_rank, num_dims, coord);


   int North,South,East,West,NW,NE,SW,SE;

   //North Neighbor
   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1];
   if(neigh[0] < 0){
      North = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &North);
   }

   //South Neighbor
   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1];
   if(neigh[0] >= process_row){
      South = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &South);
   }


   //East Neighbor
   neigh[0] = coord[0];
   neigh[1] = coord[1] + 1;
   if(neigh[1] >= process_col){
      East = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &East);
   }

   //West Neighbor
   neigh[0] = coord[0];
   neigh[1] = coord[1] - 1;
   if(neigh[1] < 0){
      West = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &West);
   }

   //North East Neighbor
   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1] + 1;
   if(neigh[0] < 0){
      NE = MPI_PROC_NULL;
   }
   else if(neigh[1] >= process_col){
      NE = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &NE);
   }

   //North West Neighbor
   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1] - 1;
   if(neigh[0] < 0){
      NW = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0){
      NW = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &NW);
   }

   //South East Neighbor
   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1] + 1;
   if(neigh[0] >= process_row){
      SE = MPI_PROC_NULL;
   }
   else if(neigh[1] >= process_col){
      SE = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &SE);
   }

   //South West Neighbor
   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1] - 1;
   if(neigh[0] >= process_row){
      SW = MPI_PROC_NULL;
   }
   else if(neigh[1] < 0){
      SW = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &SW);
   }
   



   MPI_Aint lb, extent; 
   MPI_Datatype etype, filetype, contig; 
   MPI_Offset disp;
   //basic unit of data access
   etype = MPI_UNSIGNED_CHAR;
   lb = 0; 
   MPI_Type_contiguous(cols_per_block*num_elements, MPI_UNSIGNED_CHAR, &contig);
   //extent = number of bytes to read/write + number of bytes to skip for next column
   extent = (int)processors_sqrt*cols_per_block*num_elements* sizeof(unsigned char);
   //filetype = specifies which portion of the file is visible to the process
   MPI_Type_create_resized(contig, lb, extent, &filetype);
   MPI_Type_commit(&filetype);
   //disp = number of bytes to be skipped from the start of the file 
   disp = coord[0] *(int)processors_sqrt*cols_per_block *rows_per_block*num_elements + 
      coord[1] * cols_per_block*num_elements;

   
   //Parallel I/0 read image
   
   MPI_File fh;
   // Open input image file 
   MPI_File_open(comm, filename,
      MPI_MODE_RDWR, MPI_INFO_NULL, &fh); 
   // Set view for each process
   MPI_File_set_view(fh, disp, etype, filetype, "native",
         MPI_INFO_NULL);
   //Read the bytes
   MPI_File_read(fh,Table[0],rows_per_block*cols_per_block*num_elements,
      MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
   //Close output file
   MPI_File_close(&fh);


   //Define cols datatype
   MPI_Datatype cols_type,column_type;
   MPI_Type_vector(rows_per_block,num_elements,cols_per_block*num_elements,MPI_UNSIGNED_CHAR,
      &cols_type);
   MPI_Type_commit(&cols_type);
   

   
   int tag = 11;
   //Create messages
   MPI_Request send_request[8];
   //rows
   MPI_Send_init(&(Table[0][0]),cols_per_block*num_elements, MPI_UNSIGNED_CHAR, North, 
      1, comm, &send_request[0]);
   MPI_Send_init(&(Table[rows_per_block-1][0]),cols_per_block*num_elements,
    MPI_UNSIGNED_CHAR,South,2, comm, &send_request[1]);
   //cols
   MPI_Send_init(&(Table[0][0]),1,cols_type,West,3,comm,&send_request[2]);
   MPI_Send_init(&(Table[0][(cols_per_block-1)*num_elements]),1,cols_type,East,4,comm,
      &send_request[3]);
   
   //corners
   MPI_Send_init(&(Table[0][0]),num_elements, MPI_UNSIGNED_CHAR,NW, 5, comm, 
      &send_request[4]);
   MPI_Send_init(&(Table[0][(cols_per_block-1)*num_elements]),num_elements,MPI_UNSIGNED_CHAR,NE,6,
      comm,&send_request[5]);
   MPI_Send_init(&(Table[rows_per_block-1][0]),num_elements,MPI_UNSIGNED_CHAR,SW,7,
      comm,&send_request[6]);
   MPI_Send_init(&(Table[rows_per_block-1][(cols_per_block-1)*num_elements]),num_elements,
      MPI_UNSIGNED_CHAR,SE,8, comm, &send_request[7]);



   MPI_Request receive_request[8];
   //rows
   MPI_Recv_init(halo_p->North,cols_per_block*num_elements, MPI_UNSIGNED_CHAR, North, 2, comm, 
      &receive_request[0]);
   MPI_Recv_init(halo_p->South, cols_per_block*num_elements, MPI_UNSIGNED_CHAR, South, 1, comm, 
      &receive_request[1]);
   
   //cols
   MPI_Recv_init(halo_p->West,rows_per_block*num_elements, MPI_UNSIGNED_CHAR, West, 4, comm, 
      &receive_request[2]);
   MPI_Recv_init(halo_p->East,rows_per_block*num_elements, MPI_UNSIGNED_CHAR, East, 3, comm, 
      &receive_request[3]);
   
   //corners
   MPI_Recv_init(halo_p->North_West,num_elements,MPI_UNSIGNED_CHAR,NW,8,comm,&receive_request[4]);
   MPI_Recv_init(halo_p->North_East,num_elements,MPI_UNSIGNED_CHAR,NE,7,comm,&receive_request[5]);
   MPI_Recv_init(halo_p->South_West,num_elements,MPI_UNSIGNED_CHAR,SW,6,comm,&receive_request[6]);
   MPI_Recv_init(halo_p->South_East,num_elements,MPI_UNSIGNED_CHAR,SE,5,comm,&receive_request[7]);
   

   MPI_Status status[8];
   //4x ISend
   MPI_Startall(8, send_request);
   //4x IRecv
   MPI_Startall(8, receive_request);
   MPI_Waitall(8, receive_request,status);
   
   /*
   int loop = 0;
   int changes = 0;
   int sum_changes;
   double start, finish;
   start = MPI_Wtime();
   while(loop < max_loops){
      loop++;

      //8x ISend
      MPI_Startall(8, send_request);
      //8x IRecv
      MPI_Startall(8, receive_request);
      
      //Do for our table
      for(int i=2;i<rows_per_block;i++){
         for(int j=2*num_elements;j<cols_per_block*num_elements;j++){
            if(convolution(Table, Final,i,j,h,num_elements) && !changes){
               changes++;
            }
         }
      }
      MPI_Waitall(8, receive_request,status);

      //do the job for receive
      //First row
      for(int j=num_elements;j<(cols_per_block+1)*num_elements;j++){
         if(convolution(Table, Final,1,j,h,num_elements) && !changes){
            changes++;
         }
      }
      //Last row
      for(int j=num_elements;j<(cols_per_block+1)*num_elements;j++){
         if(convolution(Table, Final,rows_per_block,j,h,num_elements) && !changes){
            changes++;
         }
      }
      //First col and last col for each middle row
      for(int i=2;i<rows_per_block;i++){
         for(int k=0;k<num_elements;k++){
            if(convolution(Table,Final,i,num_elements+k,h,num_elements) && !changes){
               changes++;
            }
            if(convolution(Table,Final,i,cols_per_block*num_elements+k,h,num_elements) && !changes){
               changes++;
            }
         }
      }

      unsigned char ** temp;
      temp = Table;
      Table = Final;
      Final = temp;
      MPI_Waitall(8,send_request,status);
      //Reduce all changes
      if(loop % 10 == 0){
         MPI_Allreduce(&changes,&sum_changes, 1, MPI_UNSIGNED_CHAR, MPI_SUM, comm);
         changes = 0;
         if(sum_changes == 0){
            if(my_rank == 0)
               printf("No changes in loop %d!\n",loop);
            break;
         }
      }
   }
   finish = MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);
   double local_elapsed, elapsed;
   local_elapsed = finish - start;
   MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   if(my_rank == 0){
      printf("Time elapsed: %f seconds\n", elapsed);
   }
   */

   // if(my_rank == 0){
   //    printf("North: ");
   //    for(int i =0;i<cols_per_block*num_elements;i++){
   //       printf("%c ",halo_p->North[i]);
   //    }
   //    printf("\nSouth: ");
   //    for(int i =0;i<cols_per_block*num_elements;i++){
   //       printf("%c ",halo_p->South[i]);
   //    }
   //    printf("\nWest: ");
   //    for(int i =0;i<rows_per_block*num_elements;i++){
   //       printf("%c ",halo_p->West[i]);
   //    }
   //    printf("\nEast: ");
   //    for(int i =0;i<rows_per_block*num_elements;i++){
   //       printf("%c ",halo_p->East[i]);
   //    }
   // }
   if (my_rank==0)
   {
      for (int i=0;i<rows_per_block;i++)
         for (int j=0;j<cols_per_block;j++)
         {
            if (j%cols_per_block==0)
               printf("\n");
            printf("%c", Table[i][j]);

         }
   }
   /* Shut down MPI */
   MPI_Barrier(MPI_COMM_WORLD);


   //MPI Parallel I/O Write the final image
   
   MPI_File fw; 
   // Open output image file
   MPI_File_open(comm, output,
      MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &fw);
   // Set view for each process
   MPI_File_set_view(fw, disp, etype, filetype, "native",
         MPI_INFO_NULL);
   //Write the bytes
   MPI_File_write(fw,Table[0],rows_per_block*cols_per_block*num_elements,
      MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
   //Close output file
   MPI_File_close(&fw); 

   Delete_Halo(halo_p);
   free(halo_p);
   free(Table[0]);
   free(Final[0]);
   free(Table);
   free(Final);
   
   MPI_Finalize(); 

   return 0;
}
