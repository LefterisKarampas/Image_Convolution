#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */ 
#include <stdlib.h>
#include <time.h>
#include <math.h>



int main(int argc,char **argv) {
   srand(time(NULL));
   int comm_sz;          
   int my_rank;
   int i;

   
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
   int num_elements = 1;
   i=1;
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
      else if ( !strcmp(argv[i], "-rgb") )
      {
         num_elements = 3;
         if(my_rank == 0){
            printf("num_elements %d\n",num_elements);
         }
         i--;
      }
      else{
         MPI_Finalize();
         return 0;
      }
      i+=2;
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
   


   unsigned char ** Table = (unsigned char **)malloc(sizeof(unsigned char *)*(rows_per_block+2));
   Table[0] = (unsigned char *)malloc(sizeof(unsigned char)*(rows_per_block+2)*(cols_per_block+2)*num_elements); 
   unsigned char ** Final = (unsigned char **)malloc(sizeof(unsigned char *)*(rows_per_block+2));
   Final[0] = (unsigned char *)malloc(sizeof(unsigned char)*(rows_per_block+2)*(cols_per_block+2)*num_elements);
   for(int i=0;i<rows_per_block+2;i++){
      Table[i]=&(Table[0][i*(cols_per_block+2)*num_elements]);
      Final[i]=&(Final[0][i*(cols_per_block+2)*num_elements]);
      if((i == 0) || (i == (rows_per_block+1))){
         for(int j=0;j<(cols_per_block+2)*num_elements;j++){
            Table[i][j] = 255;
         }
         continue;
      }
      for(int j=0;j<(cols_per_block+2)*num_elements;j++){
         if(j<num_elements){
            Table[i][j] = 255;
         }
         else if(j >=  (cols_per_block+1)*num_elements){
            Table[i][j] = 255;
         }
         else{
            Table[i][j] = my_rank;//rand() % 256; //my_rank;
         }
      }
   }
   

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
   neigh[0] = coord[0] - 1;
   neigh[1] = coord[1];
   if(neigh[0] < 0){
      North = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &North);
   }

   neigh[0] = coord[0] + 1;
   neigh[1] = coord[1];
   if(neigh[0] >= process_row){
      South = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &South);
   }

   neigh[0] = coord[0];
   neigh[1] = coord[1] + 1;
   if(neigh[1] >= process_col){
      East = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &East);
   }

   neigh[0] = coord[0];
   neigh[1] = coord[1] - 1;
   if(neigh[1] < 0){
      West = MPI_PROC_NULL;
   }
   else{
      MPI_Cart_rank(comm, neigh, &West);
   }

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
   etype = MPI_UNSIGNED_CHAR;
   lb = 0; extent = cols_per_block* sizeof(unsigned char);
   MPI_Type_contiguous(cols_per_block, MPI_UNSIGNED_CHAR, &contig);
   MPI_Type_create_resized(contig, lb, extent, &filetype);
   MPI_Type_commit(&filetype);
   disp = coord[0] * cols_per_block * processors_sqrt*num_elements + coord[1] * cols_per_block*num_elements;
   MPI_File fh,fw; 
   MPI_File_open(comm, "my_image.raw",
      MPI_MODE_RDWR, MPI_INFO_NULL, &fh); 
   MPI_File_set_view(fh, disp, etype, filetype, "native",
         MPI_INFO_NULL);
   unsigned char buf[4];
   buf[0] = -1;
   buf[1] = -2;
   MPI_File_read(fh,buf,4,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
   printf("Rank %d: %d %d %c %c\n",my_rank,(int)buf[0],(int)buf[1],buf[2],buf[3]);
   
    MPI_File_open(comm, "my_image_write.raw",
      MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &fw); 
    MPI_File_set_view(fw, disp, etype, filetype, "native",
         MPI_INFO_NULL);
   MPI_File_write(fw,buf,4,MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
   MPI_File_close(&fw); 
   MPI_File_close(&fh); 

   MPI_Barrier(MPI_COMM_WORLD);
  
   free(Table[0]);
   free(Final[0]);
   free(Table);
   free(Final);
   MPI_Finalize(); 

   return 0;
}  /* main */
