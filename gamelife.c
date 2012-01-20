/*
 *  Libraries
 */
#include <stdio.h>  /* Standard input-output                        */
#include <string.h> /* String handling                              */
#include <unistd.h> /* Miscellaneous constants, types and functions */
#include <stdlib.h> /* Standard general utilities                   */
#include <getopt.h> /* Parsing command-line options                 */
#include <time.h>   /* Get and manipulate date and time             */
#include <omp.h>    /* Shared memory multiprocessing programming    */

/*
 *  Macros
 */
#define PACKAGE      "gamelife"     /* Program name                     */
#define DEFITER      100            /* Default number of iterations     */
#define DEFROWS      9              /* Default number of rows           */
#define DEFCOLUMNS   36             /* Default number of columns        */
#define DEFINFILE    "example1.txt" /* Default input file               */
#define DEFOUTFILE   "out.txt"      /* Default output file              */
#define DEFANIMATION 0              /* Default animation state          */
#define DEFPAUSE     100            /* Default pause time               */
#define DEFMETHOD    0              /* Default method                   */
#define DEFPRINT     0              /* Default print state              */

/*
 *  Main parameters struct (go to help() implementation, line number 312)
 */
struct args_t
{
    int        iterations;   /* -i option   */
    int        rows;         /* -r option   */
    int        columns;      /* -c option   */
    char       *inFileName;  /* -f option   */
    const char *outFileName; /* -o option   */
    int        animation;    /* -a option   */
    int        pause;        /* -a argument */
    int        method;       /* -m option   */
    int        print;        /* -p option   */
};

/*
 *  Procedures and functions declaration
 */
void  gamelife(const struct args_t args);

void  initialize_world(int          **world,
                       const struct args_t args);

int   alive_neighbours(int **world,
                       int rows,
                       int columns,
                       int row,
                       int column);

void  sequential(int **world,
                 int **nextworld,
                 int rows,
                 int columns);

void  openmp(int **world,
             int **nextworld,
             int rows,
             int columns);

void  print_world(int  **world,
                  FILE *outFile,
                  int  iteration,
                  const struct args_t args);

void  help(int exitval);

int** get_memory(int rows,
                 int columns);

void  free_memory(int **matrix,
                  int rows);

/*
 * Main program
 */
int main(int argc, char* argv[])
{
  /* Variables declaration */
  struct args_t args = {.iterations = DEFITER,
                        .rows = DEFROWS,
                        .columns = DEFCOLUMNS,
                        .inFileName = DEFINFILE,
                        .outFileName = DEFOUTFILE,
                        .animation = DEFANIMATION,
                        .pause = DEFPAUSE,
                        .method = DEFMETHOD,
                        .print = DEFPRINT};
  int option;

  /* Parse command-line options */
  while((option = getopt(argc, argv, "hi:f:r:c:o:a:m:p")) != -1)
  {
    switch(option)
    {
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
        args.pause = atoi(optarg);
        break;
      case 'm':
        args.method = atoi(optarg);
        break;
      case 'p':
        args.print = 1;
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

  /* Cellular automata function*/
  gamelife(args);

  /* Main program return */
  return EXIT_SUCCESS;
}


void gamelife(const struct args_t args)
{
  int  iter;
  int  row;
  int  column;
  int  neighbours;
  int  **world;
  int  **nextworld;
  int  **tmpworld;
  FILE *outFile; 

  /* Open output file */
  if(args.print)
  {
    outFile = fopen(args.outFileName, "w");
    if(outFile == NULL)
    {
      fprintf(stderr, "%s: Error - Cannot open or create %s\n\n", PACKAGE, args.outFileName);
      help(EXIT_FAILURE);
    }
  }

  /* Allocate memory */
  world = get_memory(args.rows, args.columns);
  nextworld = get_memory(args.rows, args.columns);

  /* Initialization of first iteration */
  initialize_world(world, args); 

  if(args.print)
  {
    /* Print first iteration */
    print_world(world, outFile, -1, args);
  }

  for(iter=0; iter<args.iterations; iter++)
  {
    switch(args.method) 
    {
      case 0:
        sequential(world, nextworld, args.rows, args.columns);
        break;
      case 1:
        sequential(world, nextworld, args.rows, args.columns);
        break;
      case 2:
        //mpi(world, nextworld, args.rows, args.columns);
        break;
    }

    /* Actualize world */
    tmpworld = world;
    world = nextworld;
    nextworld = tmpworld;

    if(args.print)
    {
      /* Print iteration */
      print_world(world, outFile, iter, args);
    }
  }

  printf("\nEnd of %s\n", PACKAGE);

  /* Release memory */
  free_memory(nextworld, args.rows);
  free_memory(world, args.rows);

  if(args.print)
  {
    /* Close output file */
    fclose(outFile);
  }
}

/*
 *  World initialization procedure
 */
void initialize_world(int **world, const struct args_t args)
{
  int  i;
  int  j;
  char c, trash[100];
  FILE *inFile; 

  inFile= fopen(args.inFileName, "r");
  
  if(inFile == NULL)
  {
    fprintf(stderr, "%s: Error - File %s does not exist\n\n", PACKAGE, args.inFileName);
    help(EXIT_FAILURE);
  }
  
  i = 0;
  j = 0;
  while(!feof(inFile))
  {
    c = fgetc(inFile);
    if(c == '.') world[i][j] = 0;
    if(c == '*') world[i][j] = 1;
    j++;

    if(c == '\n' || j == args.columns)
    {
      if(j == args.columns) fgets(trash, 100, inFile);
      i++;
      j = 0;
    }

    if(i == args.rows) break; 
  }
  fclose(inFile);
}

/*
 *  Function to evaluate the number of neighbours
 */
int alive_neighbours(int **world, int rows, int columns, int row, int column)
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


/*
 *  Cellular automata sequential implementation procedure
 */
void sequential(int **world, int **nextworld, int rows, int columns)
{
  int row;
  int column;
  int neighbours;

  for(row=0; row<rows; row++)
  {
    for(column=0; column<columns; column++)
    {
      neighbours = alive_neighbours(world, rows, columns, row, column);
      if(world[row][column] == 0) nextworld[row][column] = neighbours == 3 ? 1 : 0;
      if(world[row][column] == 1) nextworld[row][column] = (neighbours == 2 || neighbours == 3) ? 1 : 0;
    }
  }
}

/*
 *  Cellular automata OpenMP implementation procedure
 */
void openmp(int **world, int **nextworld, int rows, int columns)
{
  int row;
  int column;
  int neighbours;
  
  #pragma omp parallel for private(column, neighbours)
  for(row=0; row<rows; row++)
  {
    for(column=0; column<columns; column++)
    {
      neighbours = alive_neighbours(world, rows, columns, row, column);
      if(world[row][column] == 0) nextworld[row][column] = neighbours == 3 ? 1 : 0;
      if(world[row][column] == 1) nextworld[row][column] = (neighbours == 2 || neighbours == 3) ? 1 : 0;
    }
  }
}

/*
 *  Procedure to print world content
 */
void print_world(int **world, FILE *outFile, int iteration, const struct args_t args)
{
  int  row;
  int  column;
  char c;
  char spinner[4]={'/','-','\\','|'};

  if(args.animation)
  {
    /* With animation */
    /* Restore cursor after the first iteration */
    if(iteration >= 0) printf("\033[%dA", args.rows);

    for(row=0; row<args.rows; row++)
    {
      for(column=0; column<args.columns; column++)
      {
        c = world[row][column] == 1 ? '*' : '.';
        fputc(c, outFile);
        printf("%c", c); 
      }
      fputc('\n', outFile);
      printf("\033[K\n"); // clear and next line
    }
    usleep(args.pause*1000);
  }
  else
  {
    /* Without animation */
    /* Restore cursor after the first iteration */
    printf("%c\b", spinner[(iteration+1)%sizeof(spinner)]);
    fflush(stdout);
    for(row=0; row<args.rows; row++)
    {
      for(column=0; column<args.columns; column++)
      {
        c = world[row][column] == 1 ? '*' : '.';
        fputc(c, outFile);
      }
      fputc('\n', outFile);
    }
  }

  fputc('\n', outFile);
}

/*
 *  Procedure to print help information
 */
void help(int exitval) {
  if(exitval){
    printf("%s, show working example\n", PACKAGE);
    printf("%s [-h] [-a NUM] [-f FILE] [-r NUM] [-c NUM] [-o FILE] [-m METHOD] [-p]\n\n", PACKAGE);
  }
  else {
    printf("NAME\n");
    printf("\t%s -- The game of life, introduction to parallel computing \n\n", PACKAGE);
    
    printf("SYNOPSIS\n");
    printf("\t%s [options]\n\n", PACKAGE);
    
    printf("DESCRIPTON\n");
    printf("\n\n");

    printf("OPTIONS\n");
    printf("\t%s [-h] [-a NUM] [-f FILE] [-r NUM] [-c NUM] [-o FILE] [-m METHOD] [-p]\n\n", PACKAGE);
    printf("\t-h\n\t\tprint this help and exit\n\n");
    printf("\t-a pause \n\t\tshow animation with a pause time (in miliseconds) between frames\n\n");
    printf("\t-f infile\n\t\tset intput file\n\n");
    printf("\t-r numrows\n\t\tnumber of rows\n\n");
    printf("\t-c numcolums\n\t\tnumber of columns\n\n");
    printf("\t-o outfile\n\t\tset output file\n\n");
    printf("\t-m method\n\t\tprocedure for compute next state (0=sequential, 1=OpenMP, 2=MPI)\n\n");
    printf("\t-p\n\t\tset output file\n\n");
  }
  
  exit(exitval);
}

/*
 *  Function for allocate blocks of memory from the system function
 */
int** get_memory(int rows, int columns)
{
  int i;
  int **matrix;
  
  matrix = (int **) malloc(rows * sizeof(int *));

  if(matrix == (int **) NULL)
  {
    fprintf(stderr, "%s: Error - Memory error\n\n", PACKAGE);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<rows; i++)
  {
    matrix[i] = (int *) malloc(columns * sizeof(int));

    if(matrix[i] == (int *) NULL)
    {
      fprintf(stderr, "%s: Error - Memory error\n\n", PACKAGE);
      exit(EXIT_FAILURE);
    }
    else memset(matrix[i], 0, columns);
  }

  return matrix;
}

/*
 *  Function for release blocks of memory back to the system
 */
void free_memory(int **matrix, int rows)
{
  int i;

  for(i=0; i<rows; i++) free(matrix[i]);

  free(matrix);
}
