#include <stdio.h>              /* Entrada/Salida */
#include <string.h>             /* Operaciones con cadenas */
#include <unistd.h>             /* Misceláneo (constantes, tipos y funciones) */
#include <stdlib.h>             /* Manejo de memoria dinámica, procesos, etc */
#include <getopt.h>             /* Manejo de argumentos */
#include <time.h>
#include <omp.h>

#define PACKAGE "gamelife"
#define DEFITER 100
#define DEFROWS 9
#define DEFCOLUMNS 36
#define DEFINFILE "example1.txt"
#define DEFOUTFILE "out.txt"



struct args_t {
    int iterations;             /* -i option */
    int rows;                   /* -r option */
    int columns;                /* -c option */
    char *inFileName;           /* -f option */
    const char *outFileName;    /* -o option */
    int animation;              /* -a option */
    int method;                 /* -m option */
};


void gamelife(const struct args_t args);
void initialize_world(int **world, const struct args_t args);
int alive_neighbors(int **world, int rows, int columns, int row, int column);
void print_world(int **world, int rows, int columns, FILE *outFile, int animation, int iteration);
void sequential(int **world, int **nextworld, int rows, int columns);
void openmp(int **world, int **nextworld, int rows, int columns);
void help(int exitval);
int** get_memory(int rows, int columns);
void free_memory(int **matrix, int rows);



int main(int argc, char* argv[]) {
  struct args_t args = {.iterations = DEFITER,.rows = DEFROWS, .columns = DEFCOLUMNS, .inFileName = DEFINFILE, .outFileName = DEFOUTFILE, .animation = 0, .method = 0};
  int option;

  while((option = getopt(argc, argv, "hi:f:r:c:o:am:")) != -1) {
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
      case ':':
        fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE, optopt);
        help(EXIT_FAILURE);
        break;
      case '?':
        fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
        help(1);
    }
  }

  /* Print all remaining options */
  for(; optind < argc; optind++)
    printf("Non-option argument: %s\n", argv[optind]);

  gamelife(args);

  return EXIT_SUCCESS;
}


void sequential(int **world, int **nextworld, int rows, int columns) {
  int row, column, neighbors;

  for(row=0; row<rows; row++) {
    for(column=0; column<columns; column++) {
      neighbors = alive_neighbors(world, rows, columns, row, column);
      if(world[row][column] == 0) 
        nextworld[row][column] = neighbors == 3 ? 1 : 0;

      if(world[row][column] == 1)
        nextworld[row][column] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    }
  }
}


void openmp(int **world, int **nextworld, int rows, int columns) {
  int row, column, neighbors;
  
  #pragma omp parallel for private(column, neighbors) if(rows >= columns)
  for(row=0; row<rows; row++) {
    #pragma omp parallel for private(row, neighbors) if(columns > rows)
    for(column=0; column<columns; column++) {
      neighbors = alive_neighbors(world, rows, columns, row, column);
      if(world[row][column] == 0) 
        nextworld[row][column] = neighbors == 3 ? 1 : 0;

      if(world[row][column] == 1)
        nextworld[row][column] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    }
  }
}


void gamelife(const struct args_t args) {
  int iter, row, column, neighbors;
  int **world, **nextworld, **tmpworld;
  FILE *outFile; 

  outFile = fopen(args.outFileName, "w");
  
  if(outFile == NULL) {
    fprintf(stderr, "%s: Error - Cannot open or create %s\n\n", PACKAGE, args.outFileName);
    help(EXIT_FAILURE);
  }

  world = get_memory(args.rows, args.columns);
  nextworld = get_memory(args.rows, args.columns);

  initialize_world(world, args); 

  print_world(world, args.rows, args.columns, outFile, args.animation, -1);

  for(iter=0; iter<args.iterations; iter++) {
    if(args.method == 0)
      sequential(world, nextworld, args.rows, args.columns);
    
    if(args.method ==1)
      openmp(world, nextworld, args.rows, args.columns);

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


void initialize_world(int **world, const struct args_t args) {
  int i, j;
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
    
    if(c == '.')
      world[i][j] = 0;
    
    if(c == '*')
      world[i][j] = 1;

    j++;

    if(c == '\n' || j == args.columns) {
      if(j == args.columns)
        fgets(trash, 100, inFile);
      i++;
      j = 0;
    }

    if(i == args.rows)
      break; 
  }
  fclose(inFile);
}


int alive_neighbors(int **world, int rows, int columns, int row, int column) {
  int rowLeft, rowRight, columnUp, columnBottom;
  rowLeft = (row - 1 + rows)%rows;
  rowRight = (row + 1 + rows)%rows;
  columnUp = (column - 1 + columns)%columns;
  columnBottom = (column + 1 + columns)%columns;
  
  return  world[rowLeft][columnUp] + 
          world[rowLeft][column] + 
          world[rowLeft][columnBottom] + 
          world[row][columnUp] + 
          world[row][columnBottom] + 
          world[rowRight][columnUp] + 
          world[rowRight][column] + 
          world[rowRight][columnBottom]; 
}


void print_world(int **world, int rows, int columns, FILE *outFile, int animation, int iteration) {
  int row, column;
  char c;
  char spinner[4]={'/','-','\\','|'};
  
  if(animation && iteration >= 0 )
    printf("\033[%dA", rows); // restore cursor after the first iteration

  if(!animation) {
    printf("%c\b", spinner[(iteration+1)%sizeof(spinner)]);
    fflush(stdout);
  }

  for(row=0; row<rows; row++) {
    for(column=0; column<columns; column++) {
      c = world[row][column] == 1 ? '*' : '.';
      fputc(c, outFile);
      if(animation)
        printf("%c", c); 
    }
    fputc('\n', outFile);
    if(animation)
       printf("\033[K\n"); // clear and next line
  }
      
  if(animation)
    usleep(200000);
    
  fputc('\n', outFile);
}


void help(int exitval) {
  printf("NAME\n");
  printf("%s -- \n", PACKAGE); 
  printf("%s [-h] [-a] [-f FILE] [-r NUM] [-c NUM] [-o FILE] [-m METHOD]\n\n", PACKAGE);
  printf("  -h               print this help and exit\n");
  printf("  -a               show animation\n\n");
  printf("  -f FILE          set intput file\n");
  printf("  -r NUM           number of rows\n");
  printf("  -c NUM           number of columns\n");
  printf("  -o FILE          set output file\n\n");
  
  exit(exitval);
}


int** get_memory(int rows, int columns) {
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
    else
      memset(matrix[i], 0, columns);
  }

  return matrix;
}


void free_memory(int **matrix, int rows){
  int i;

  for(i=0; i<rows; i++)
    free(matrix[i]);

  free(matrix);
}
