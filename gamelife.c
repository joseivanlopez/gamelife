/*
Opciones por añadir:
  -cambiar tiempo de pausa entre iteracion e iteración en la animación.
  -añadir opción para que guarde en fichero o no.
  -añadir cálculo de tiempos de ejecucion
*/

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
#define DEFMETHOD    0              /* Default method                   */

/*
 *  Main parameters struct (go to help() implementation, line number 266)
 */
struct args_t
{
    int        iterations;   /* -i option */
    int        rows;         /* -r option */
    int        columns;      /* -c option */
    char       *inFileName;  /* -f option */
    const char *outFileName; /* -o option */
    int        animation;    /* -a option */
    int        method;       /* -m option */
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

void  print_world(int  **world,
                  int  rows,
                  int  columns,
                  FILE *outFile,
                  int  animation,
                  int  iteration);

void  sequential(int **world,
                 int **nextworld,
                 int rows,
                 int columns);

void  openmp(int **world,
             int **nextworld,
             int rows,
             int columns);

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
                        .method = DEFMETHOD};
  int option;

  /* Parse command-line options */
  while((option = getopt(argc, argv, "hi:f:r:c:o:am:")) != -1)
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
        break;
      case 'm':
        args.method = atoi(optarg);
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
  
  #pragma omp parallel for private(column, neighbours)// if(rows >= columns)
  for(row=0; row<rows; row++)
  {
    //#pragma omp parallel for private(row, neighbours)// if(columns > rows)
    for(column=0; column<columns; column++)
    {
      neighbours = alive_neighbours(world, rows, columns, row, column);
      if(world[row][column] == 0) nextworld[row][column] = neighbours == 3 ? 1 : 0;
      if(world[row][column] == 1) nextworld[row][column] = (neighbours == 2 || neighbours == 3) ? 1 : 0;
    }
  }
}

/*
 *  Cellular automata procedure calling
 */
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

  outFile = fopen(args.outFileName, "w");
  
  if(outFile == NULL)
  {
    fprintf(stderr, "%s: Error - Cannot open or create %s\n\n", PACKAGE, args.outFileName);
    help(EXIT_FAILURE);
  }

  world = get_memory(args.rows, args.columns);
  nextworld = get_memory(args.rows, args.columns);

  initialize_world(world, args); 

  print_world(world, args.rows, args.columns, outFile, args.animation, -1);

  for(iter=0; iter<args.iterations; iter++)
  {
    if(args.method == 0) sequential(world, nextworld, args.rows, args.columns);
    if(args.method == 1) openmp(world, nextworld, args.rows, args.columns);

    tmpworld = world;
    world = nextworld;
    nextworld = tmpworld;
    
    print_world(world, args.rows, args.columns, outFile, args.animation, iter);
  }

  printf("\nEnd of %s\n", PACKAGE);

  free_memory(nextworld, args.rows);
  free_memory(world, args.rows);
  fclose(outFile);
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
 *  Procedure to print world content
 */
void print_world(int **world, int rows, int columns, FILE *outFile, int animation, int iteration)
{
  int  row;
  int  column;
  char c;
  char spinner[4]={'/','-','\\','|'};

  /* Restore cursor after the first iteration */
  if(animation && iteration >= 0) printf("\033[%dA", rows);

  if(!animation)
  {
    printf("%c\b", spinner[(iteration+1)%sizeof(spinner)]);
    fflush(stdout);
  }

  for(row=0; row<rows; row++)
  {
    for(column=0; column<columns; column++)
    {
      c = world[row][column] == 1 ? '*' : '.';
      fputc(c, outFile);
      if(animation) printf("%c", c); 
    }
    fputc('\n', outFile);
    if(animation) printf("\033[K\n"); // clear and next line
  }
      
  if(animation) usleep(200000); //0-1000000 microseconds
    
  fputc('\n', outFile);
}

/*
 *  Procedure to print help information
 */
void help(int exitval)
{
  printf("NAME\n");
  printf("%s -- \n", PACKAGE); 
  printf("%s [-h] [-a] [-f FILE] [-r NUM] [-c NUM] [-o FILE] [-m METHOD]\n\n", PACKAGE);
  printf("  -h               print this help and exit\n");
  printf("  -a               show animation\n\n");
  printf("  -f FILE          set intput file\n");
  printf("  -r NUM           number of rows\n");
  printf("  -c NUM           number of columns\n");
  printf("  -o FILE          set output file\n\n");
  printf("  -m NUM           process method (0 = sequential, 1 = OpenMP, 2 = MPI)\n\n");
  
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
