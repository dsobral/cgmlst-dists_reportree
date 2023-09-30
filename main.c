#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/types.h>
#include <inttypes.h>
#include <pthread.h>
#include <math.h>

#define VERSION "0.4.0"
#define EXENAME "cgmlst-dists"
#define GITHUB_URL "https://github.com/genpat-it/cgmlst-dists-64"
//#define DEBUG

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
	 
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
	 
const int32_t MAX_LINE = 1E5;
const int32_t MAX_ASM  = 1E5;
const char* DELIMS = "\n\r\t";
const uint8_t IGNORE_ALLELE = 0;
const char REPLACE_CHAR = ' ';


struct targs {
  int32_t quiet;
  int64_t thread;
  int64_t t1;
  int64_t t2;	
  uint32_t** call;
  uint32_t* dist;
  //int64_t j;
  int64_t nrow;
  size_t ncol;
  uint32_t maxdiff;
};


//------------------------------------------------------------------------
void show_help(uint32_t retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);

  static const char str[] = {
      "SYNOPSIS\n  Pairwise CG-MLST distance matrix from allele call tables\n"
      "USAGE\n  %s [options] chewbbaca.tab > distances.tsv\n"
      "OPTIONS\n"
      "  -h\tShow this help\n"
      "  -v\tPrint version and exit\n"
      "  -q\tQuiet mode; do not print progress information\n"
      "  -c\tUse comma instead of tab in output\n"
      "  -m N\tOutput: 1=lower-tri 2=upper-tri 3=full [3]\n"
      "  -x N\tStop calculating beyond this distance [9999]\n"
      "  -t N\tNumber of threads to use [1]\n"
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
uint32_t distance(const uint32_t* restrict a, const uint32_t* restrict b, size_t len, uint32_t maxdiff)
{
  uint32_t diff=0;
  for (size_t i=0; i < len; i++) {
    if (a[i] != b[i] && a[i] != IGNORE_ALLELE && b[i] != IGNORE_ALLELE) {
      diff++;
      if (diff >= maxdiff) return maxdiff;
    }
  }
  return diff;
}

void *thread_function(void *args)
{

  struct targs *values;
  values = (struct targs *)args;  
  
  int32_t quiet = values->quiet;
  
  int64_t thread = values->thread;
  
  int64_t t1 = values->t1;
  int64_t t2 = values->t2;
  
  fprintf(stderr, "Thread %ld is running from row %ld to %ld\n",thread,t1,t2);
	 
  uint32_t** call = values->call;
  uint32_t* dist = values->dist;

  int64_t nrow = values->nrow;
  size_t ncol = values->ncol;
  uint32_t maxdiff = values->maxdiff;
  
  for (int64_t j=0; j < (t2-t1); j++) {	 
    if (!quiet) fprintf(stderr, "\rThread %ld working: %.2f%%", thread, (j+1)*100.0/(t2-t1));  
	for (int64_t i=0; i < nrow; i++) {
		uint32_t d = distance(call[j], call[i], ncol, maxdiff);
		dist[j*nrow+i] = d;
	}
  }
  fprintf(stderr, "Thread %ld finished\n",values->thread);
}

//------------------------------------------------------------------------
void* calloc_safe(size_t nmemb, size_t size) 
{
  void* ptr = calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "ERROR: could not allocate %ld kb RAM\n", (nmemb*size)>>10);
    exit(EXIT_FAILURE);
  }
  return ptr;
}

//------------------------------------------------------------------------

uint32_t str_replace(char* str, char* old, char* new)
{
  size_t sl = strlen(str);
  size_t ol = strlen(old);
  size_t nl = strlen(new);
  if (ol < 1 || nl < 1 || ol != nl || sl < ol) {
    fprintf(stderr, "ERROR: str_replace(%lu,%lu,%lu)\n", sl, ol, nl);
    exit(EXIT_FAILURE);
  }

  // char *strstr(const char *haystack, const char *needle);
  char* p = NULL;
  while ( (p = strstr(str, old)) != NULL ) {
    // char *strncpy(char *dest, const char *src, size_t n);
    strncpy(p, new, nl);
  }
  
  return sl;
}


//------------------------------------------------------------------------

void cleanup_line(char* str)
{
  char* s = str;
#ifdef DEBUG
  fprintf(stderr, "BEFORE: %s", str);
#endif
  // skip over first column (ID)
  while (*s != 0 && *s != '\t') s++;

  // CHewbacca codes
  // LINF NIPH INF-nnn PLOT3 PLOT5 ASM
  // these two are special as they end in numbers
  // and don't want to confuse with INF-xxx
  str_replace(s, "PLOT3", "    0");
  str_replace(s, "PLOT5", "    0");

  // replace alpha with space so atoi() works
  while (*s++) {
    if (isalpha(*s)) {
      *s = REPLACE_CHAR;
    }
  }
#ifdef DEBUG
  fprintf(stderr, "AFTER : %s", str);
#endif
}

//------------------------------------------------------------------------
int32_t main(int argc, char* argv[])
{
  // parse command line parameters
  int32_t opt, quiet = 0, csv = 0, threads = 1, mode = 3, maxdiff = 9999;
  while ((opt = getopt(argc, argv, "hqcvm:t:x:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
      case 't': threads = atoi(optarg); break;
      case 'm': mode = atoi(optarg); break;
      case 'x': maxdiff = atoi(optarg); break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      default: show_help(EXIT_FAILURE);
    }
  }

  // require a filename argument
  if (optind >= argc) {
    show_help(EXIT_FAILURE);
    return 0;
  }
  const char* infile = argv[optind];

  // parameter sanity check
  if (threads < 1) threads = 1;

  // say hello
  if (!quiet) {
    fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);
  }

  // read file one line at a time
  FILE* in = (FILE*) fopen(infile, "r");
  if (! in) {
    fprintf(stderr, "ERROR: can not open file '%s'\n", infile);
    exit(EXIT_FAILURE);
  }

  // allocate RAM
  char* buf = (char*) calloc_safe( MAX_LINE, sizeof(char) );

  char** id  = (char**) calloc_safe( MAX_ASM, sizeof(char*) );
  uint32_t** call = (uint32_t**) calloc_safe( MAX_ASM, sizeof(uint32_t*) );

  int32_t row = -1;
  uint32_t ncol = 0;
   
  while (fgets(buf, MAX_LINE, in))
  {
    // cleanup non-numerics in NON-HEADER lines
    if (row >=0) cleanup_line(buf);

    // scan for tab separated values
    char* save;
    char* s = strtok_r(buf, DELIMS, &save);
    int32_t col = -1;
    while (s) {
      if (row >= 0) {
        if (col < 0) {
          if (strlen(s)==0) {
            fprintf(stderr, "row %d has an empty ID in first column\n", row+1);
            exit(EXIT_FAILURE);
          }
          id[row] = strdup(s);
          call[row] = (uint32_t*) calloc_safe(ncol, sizeof(uint32_t*));
        }
        else {
          // INF-xxxx are returned as -ve numbers
          call[row][col] = abs( atoi(s) );
        }
      }
      else {
        // just parsing columns on first row
      }
      col++;
      s = strtok_r(NULL, DELIMS, &save);
    }
    row++;
//    if (!quiet) fprintf(stderr, "row %d has %d cols\n", row, col) ;
    if (!quiet) fprintf(stderr, "\rLoaded row %d", row);
    if (row==0) ncol = col;
    if (row >= MAX_ASM) {
      fprintf(stderr, "Too many rows, can only handle %d\n", MAX_ASM);
      exit(EXIT_FAILURE);
    }
    
    if (col != ncol) {
      fprintf(stderr, "\nERROR: row %d had %d cols, expected %d\n", row+1, col+1, ncol);
      exit(-1);
    }
  }
  int64_t nrow = row;
  fclose(in);
  
  // what we collected
  if (!quiet) fprintf(stderr, "\rLoaded %ld samples x %d allele calls\n", nrow, ncol);  
   	
  uint32_t *dist_list[threads];
  
  pthread_t tids[threads];
  struct targs *argset[threads];
	
  int64_t interval = nrow / threads;
  int64_t interval_remainder = nrow % threads;
	
  for(int64_t t=0;t<threads;t++){
		
	struct targs *th = (struct targs *)malloc(sizeof(struct targs));
	
	th->quiet = quiet;
	th->thread = t;
	th->t1 = t*interval;
	th->t2 = (t+1)*interval;
	
	//Last thread will take whatever remains...
	if(t==(threads-1)){
		th->t2 = th->t2 + interval_remainder;
	}

	uint32_t *thread_dist = calloc_safe((th->t2 - th->t1)*nrow, sizeof(uint32_t));	
	th->dist = thread_dist;
	
	dist_list[t] = thread_dist;
	
	th->call = call;

	th->nrow = nrow;
	th->ncol = ncol;
	th->maxdiff = maxdiff;
		
	pthread_t newThread;
	pthread_create(&newThread, NULL, thread_function, (void *)th);
	tids[t] = newThread;
	argset[t] = th;
  }
	
  // wait for them to finish... 
  // hopefully they don't have to wait much after one finishes...
  for(int64_t t=0;t<threads;t++){
	pthread_join(tids[t], NULL);
  }
	
  if (!quiet) fprintf(stderr, "\nWriting distance matrix to stdout...\n");

  // separator choice
 
  // Print header row
  char sep = csv ? ',' : '\t';
  for (uint32_t j=0; j < nrow; j++) {
    if (j==0) printf(EXENAME);
    printf("%c%s", sep, id[j]);
  }
  printf("\n");

  // Print matrix
  for(int64_t t=0;t<threads;t++){  
    for (int64_t j=0; j < (argset[t]->t2 - argset[t]->t1) ; j++) {
      printf("%s", id[t*interval + j]);
      uint32_t start = (mode & 1) ?    0 : j   ;  // upper?
      uint32_t end   = (mode & 2) ? nrow : j+1 ;  // lower?
      for (int64_t i=start; i < end; i++) {
        printf("%c%d", sep, dist_list[t][j*nrow + i]);
      }
      printf("\n");
    }
  }

  // free RAM
  for(int64_t t=0;t<threads;t++){
	free(argset[t]);
	free(dist_list[t]);
  }
  for (uint32_t i=0; i < nrow; i++) {
    free( id[i] );
    free( call[i] );
  }
  free(id);
  free(call);
  free(buf);

  if (!quiet) fprintf(stderr, "\nDone.\n");

  return 0;
}

//------------------------------------------------------------------------