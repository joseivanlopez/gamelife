#include <stdio.h>  /* Standard input-output                        */
#include <string.h> /* String handling                              */
#include <unistd.h> /* Miscellaneous constants, types and functions */
#include <stdlib.h> /* Standard general utilities                   */
#include <getopt.h> /* Parsing command-line options                 */
#include <time.h>   /* Get and manipulate date and time             */
#include <omp.h>    /* Shared memory multiprocessing programming    */
#include <mpi.h>

#define PACKAGE       "gamelife"     /* Program name                     */
#define DEFITER       100            /* Default number of iterations     */
#define DEFROWS       9              /* Default number of rows           */
#define DEFCOLUMNS    36             /* Default number of columns        */
#define DEFINFILE     "example1.txt" /* Default input file               */
#define DEFOUTFILE    "out.txt"      /* Default output file              */
#define DEFMETHOD     0              /* Default method                   */
#define DEFANIMATION  0              /* Default animation state          */
#define DEFSTEP       100            /* Default pause time               */
#define DEFTIMES      0

/* Struct for program arguments */
struct args_t
{
  int iterations;     /* -i option  */
  int rows;           /* -r option  */
  int columns;        /* -c option  */
  char *inFileName;   /* -f option  */
  char *outFileName;  /* -o option  */
  int method;         /* -m option  */
  int animation;      /* -a option  */
  int step;           /* -s option  */
  int times;          /* -t option  */
};

/* Procedures and functions declaration */
void gamelife(const struct args_t args);
void sequential(int **world, int **nextworld, int rows, int columns);
void openmp(int **world, int **nextworld, int rows, int columns);
void mpi(int **world, int **nextworld, int rows, int columns, int root);
void initialize_world(int **world, const struct args_t args);
int alive_neighbors(int **world, int rows, int columns, int row, int column);
void print_world(int **world, FILE *outFile, int iteration, const struct args_t args);
void help(int exitval);
int** get_memory(int rows, int columns);
void free_memory(int **matrix, int rows);


/* Main program */
int main(int argc, char* argv[])
{
  struct args_t args = { 
    .iterations = DEFITER,
    .rows = DEFROWS,
    .columns = DEFCOLUMNS,
    .inFileName = DEFINFILE,
    .outFileName = DEFOUTFILE,
    .method = DEFMETHOD,
    .animation = DEFANIMATION,
    .step = DEFSTEP,
    .times = DEFTIMES
  };
  int option;

  /* Parse command-line options */
  while((option = getopt(argc, argv, "hi:r:c:f:o:m:as:t")) != -1) {
    switch(option) {
      case 'h':
        help(EXIT_SUCCESS);
        break;
      case 'i':
        args.iterations = atoi(optarg);
        break;
      case 'r':
        args.rows = atoi(optarg);
        break;
      case 'c':
        args.columns = atoi(optarg);
        break;
      case 'f':
        args.inFileName = optarg;
        break;
      case 'o':
        args.outFileName = optarg;
        break;
      case 'a':
        args.animation = 1;
        break;
      case 'm':
        args.method = atoi(optarg);
        break;
      case 's':
        args.step = atoi(optarg);
        break;
      case 't':
        args.times = 1;
        break;
      case ':':
        fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE, optopt);
        help(EXIT_FAILURE);
        break;
      case '?':
        fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
        help(1);
    }
  }
  for(; optind < argc; optind++) printf("Non-option argument: %s\n", argv[optind]);
  
  switch(args.method) {
    case 2: 
      MPI_Init(&argc, &argv);
      gamelife(args);
      MPI_Finalize();
      break;
    default:
      gamelife(args);
  }

  return EXIT_SUCCESS;
}


void gamelife(const struct args_t args) {
  int iter, row, column, neighbors;
  int **world, **nextworld, **tmpworld;
  FILE *outFile;
  double start, stop; 
  int root, rank;
  
  root = 0;  

  if(args.method == 2) 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  else 
    rank = 0;

  /* Open output file */
  if(!args.times && rank == root) {
    outFile = fopen(args.outFileName, "w");
    if(outFile == NULL) {
      fprintf(stderr, "%s: Error - Cannot open or create %s\n\n", PACKAGE, args.outFileName);
      help(EXIT_FAILURE);
    }
  }

  /* Allocate memory */
  world = get_memory(args.rows, args.columns);
  nextworld = get_memory(args.rows, args.columns);

  /* Initialization of first iteration */
  initialize_world(world, args); 

  /* Print first iteration */
  if(!args.times && rank == root) {
    print_world(world, outFile, -1, args);
  }

  start = omp_get_wtime();

  for(iter=0; iter<args.iterations; iter++) {
    switch(args.method) {
      case 0:
        sequential(world, nextworld, args.rows, args.columns);
        break;
      case 1:
        openmp(world, nextworld, args.rows, args.columns);
        break;
      case 2:
        mpi(world, nextworld, args.rows, args.columns, root);
        break;
    }

    /* Actualize world */
    tmpworld = world;
    world = nextworld;
    nextworld = tmpworld;

    /* Print iteration */
    if(!args.times && rank == root) {
      print_world(world, outFile, iter, args);
    }
  }

  stop = omp_get_wtime();

  if(args.times && rank == root) {
    printf("%d;%d;%d;%d;%.3f\n", args.method, args.rows, args.columns, args.iterations, stop-start);
  }
  
  if(!args.times && rank == root)  
    printf("\nEnd of %s.\n", PACKAGE);

  /* Release memory */
  free_memory(nextworld, args.rows);
  free_memory(world, args.rows);

  /* Close output file */
  if(!args.times && rank == root) {
    fclose(outFile);
  }
}


/* Cellular automata sequential implementation procedure */
void sequential(int **world, int **nextworld, int rows, int columns) {
  int row, column, neighbors;

  for(row=0; row<rows; row++) {
    for(column=0; column<columns; column++) {
      neighbors = alive_neighbors(world, rows, columns, row, column);
      if(world[row][column] == 0) nextworld[row][column] = neighbors == 3 ? 1 : 0;
      if(world[row][column] == 1) nextworld[row][column] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    }
  }
}

/* Cellular automata OpenMP implementation procedure */
void openmp(int **world, int **nextworld, int rows, int columns) {
  int row, column, neighbors;
  
  #pragma omp parallel for private(column, neighbors)
  for(row=0; row<rows; row++) {
    for(column=0; column<columns; column++) {
      neighbors = alive_neighbors(world, rows, columns, row, column);
      if(world[row][column] == 0) nextworld[row][column] = neighbors == 3 ? 1 : 0;
      if(world[row][column] == 1) nextworld[row][column] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    }
  }
}


void mpi(int **world, int **nextworld, int rows, int columns, int root) {
  int rank, prank, prevrank, nextrank, numprocs, np;
  int first, last, row, column, numrows, size, neighbors;
  int **recvbuf, **sendbuf;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  
  np = (numprocs >= rows) ? rows : numprocs; 

  /* root envía una parte de la matriz a cada proceso */
  if(rank == root) {
    numrows = rows/np;
    for(prank=0; prank<np; prank++) {
      first = prank*numrows;
      last = first + numrows;
      last = last > rows ? rows : last;
      size = last - first;
      /* Número de filas que se le va a enviar a cada proceso */
      MPI_Send(&size, 1, MPI_INT, prank, 0, MPI_COMM_WORLD);
      for(row=first; row<last; row++) {
       /* Se envía una fila */
        MPI_Send(world[row], columns, MPI_INT, prank, 0, MPI_COMM_WORLD);
      }
    }
  }  

  /* Recibe el número de filas */
  MPI_Recv(&numrows, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
  
  /* Recibe cada fila */
  recvbuf = get_memory(numrows+2, columns);
  for(row=1; row<=numrows; row++) {
    MPI_Recv(recvbuf[row], columns, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
  }
  
  /* Cada proceso envía su primera y última fila y recibe las del resto */
  prevrank = (rank-1+np)%np;
  nextrank = (rank+1)%np; 
  
  MPI_Send(recvbuf[1], columns, MPI_INT, prevrank, 0, MPI_COMM_WORLD);
  MPI_Send(recvbuf[numrows], columns, MPI_INT, nextrank, 0, MPI_COMM_WORLD);
  MPI_Recv(recvbuf[0], columns, MPI_INT, prevrank, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(recvbuf[numrows+1], columns, MPI_INT, nextrank, 0, MPI_COMM_WORLD, &status);

  /* Cada proceso realiza los cálculos con su parte de la matriz */
  sendbuf = get_memory(numrows, columns);
  for(row=1; row<=numrows; row++) {
    for(column=0; column<columns; column++) {
      neighbors = alive_neighbors(recvbuf, numrows+2, columns, row, column);
      if(recvbuf[row][column] == 0) sendbuf[row-1][column] = neighbors == 3 ? 1 : 0;
      if(recvbuf[row][column] == 1) sendbuf[row-1][column] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    }
    MPI_Send(sendbuf[row-1], columns, MPI_INT, root, 0, MPI_COMM_WORLD);
  }  

  /* root recibe cada fila modificada */ 
  if(rank == root) {
    for(prank=0; prank<np; prank++) {
      first = prank*size;
      last = first + numrows;
      last = last > rows ? rows : last;
      for(row=first; row<last; row++) {
        MPI_Recv(nextworld[row], columns, MPI_INT, prank, 0, MPI_COMM_WORLD, &status);
      }
    }
  }
      
  free_memory(recvbuf, numrows+2);
  free_memory(sendbuf, numrows);
}


/* World initialization procedure */
void initialize_world(int **world, const struct args_t args)
{
  int  i, j;
  char c, trash[100];
  FILE *inFile; 

  inFile= fopen(args.inFileName, "r");
  
  if(inFile == NULL) {
    fprintf(stderr, "%s: Error - File %s does not exist\n\n", PACKAGE, args.inFileName);
    help(EXIT_FAILURE);
  }
  
  i = 0;
  j = 0;
  while(!feof(inFile)) {
    c = fgetc(inFile);
    if(c == '.') world[i][j] = 0;
    if(c == '*') world[i][j] = 1;
    j++;

    if(c == '\n' || j == args.columns) {
      if(j == args.columns) fgets(trash, 100, inFile);
      i++;
      j = 0;
    }

    if(i == args.rows) break; 
  }
  fclose(inFile);
}


/* Function to evaluate the number of neighbors */
int alive_neighbors(int **world, int rows, int columns, int row, int column)
{
  int rowLeft      = (row - 1 + rows)%rows;
  int rowRight     = (row + 1 + rows)%rows;
  int columnUp     = (column - 1 + columns)%columns;
  int columnBottom = (column + 1 + columns)%columns;
  
  return world[rowLeft][columnUp] + 
         world[rowLeft][column] + 
         world[rowLeft][columnBottom] + 
         world[row][columnUp] + 
         world[row][columnBottom] + 
         world[rowRight][columnUp] + 
         world[rowRight][column] + 
         world[rowRight][columnBottom]; 
}


/* Procedure to print world content */
void print_world(int **world, FILE *outFile, int iteration, const struct args_t args)
{
  int row, column;
  char c;
  char spinner[4]={'/','-','\\','|'};

  /* Restore cursor after the first iteration */
  if(args.animation && iteration >= 0) printf("\033[%dA", args.rows);
  
  if(!args.animation) {
    printf("%c\b", spinner[(iteration+1)%sizeof(spinner)]);
    fflush(stdout);
  }

  for(row=0; row<args.rows; row++) {
    for(column=0; column<args.columns; column++) {
      c = world[row][column] == 1 ? '*' : '.';
      fputc(c, outFile);
      if(args.animation) printf("%c", c); 
    }
    fputc('\n', outFile);
    if(args.animation) printf("\033[K\n"); // clear and next line
  }
  fputc('\n', outFile);
  if(args.animation) usleep(args.step*1000);
}


/* Procedure to print help information */
void help(int exitval) {
  if(exitval){
    printf("%s, show working example\n", PACKAGE);
    printf("%s [-h] [-f FILE] [-o FILE] [-m METHOD] [-r NUM] [-c NUM] [-a] [-s STEP] [-t]\n\n", PACKAGE);
  }
  else {
    printf("NAME\n");
    printf("\t%s -- The game of life, introduction to parallel computing \n\n", PACKAGE);
    
    printf("SYNOPSIS\n");
    printf("\t%s [options]\n\n", PACKAGE);
    
    printf("DESCRIPTON\n");
    printf("\n\n");

    printf("OPTIONS\n");
    printf("\t%s [-h] [-f FILE] [-o FILE] [-m METHOD] [-r NUM] [-c NUM] [-a] [-s STEP] [-t]\n\n", PACKAGE);
    printf("\t-h\n\t\tprint this help and exit\n\n");
    printf("\t-f infile\n\t\tset intput file\n\n");
    printf("\t-o outfile\n\t\tset output file\n\n");
    printf("\t-m method\n\t\tprocedure for compute next state (0=sequential, 1=OpenMP, 2=MPI)\n\n");
    printf("\t-r numrows\n\t\tnumber of rows\n\n");
    printf("\t-c numcolums\n\t\tnumber of columns\n\n");
    printf("\t-a \n\t\tshow animation with a pause time (in miliseconds) between frames\n\n");
    printf("\t-s\n\t\tset step time for animation\n\n");
    printf("\t-t\n\t\ttime measure\n\n");
    
    printf("EXAMPLES\n");
    printf("\tWith Sequential or openMP mode:\n");
    printf("\t$ ./gamelife [options] -m [0|1]\n");
    printf("\t$ OMP_NUM_THREADS=N ./gamelife [options] -m 1\n\n");
    printf("\tWith MPI mode:\n");
    printf("\t$ mpirun -np N gamelife [options] -m 2\n\n");
  }
  
  exit(exitval);
}


/* Function for allocate blocks of memory from the system function */
int** get_memory(int rows, int columns)
{
  int i;
  int **matrix;
  
  matrix = (int **) malloc(rows * sizeof(int *));

  if(matrix == (int **) NULL) {
    fprintf(stderr, "%s: Error - Memory error\n\n", PACKAGE);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<rows; i++) {
    matrix[i] = (int *) malloc(columns * sizeof(int));

    if(matrix[i] == (int *) NULL) {
      fprintf(stderr, "%s: Error - Memory error\n\n", PACKAGE);
      exit(EXIT_FAILURE);
    }
    else memset(matrix[i], 0, columns);
  }

  return matrix;
}

/* Function for release blocks of memory back to the system */
void free_memory(int **matrix, int rows) {
  int i;

  for(i=0; i<rows; i++) free(matrix[i]);

  free(matrix);
}
