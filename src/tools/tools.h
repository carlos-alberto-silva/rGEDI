#define TIDY(arr) if((arr)){free((arr));(arr)=NULL;}  /*free an array*/

double *dalloc(int,char *,int);
float *falloc(uint64_t,char *,int);
char *challoc(uint64_t,char *,int);
unsigned char *uchalloc(uint64_t,char *,int);
int *ialloc(int,char *,int);
short int *shalloc(int,char *,int);
int *intSwap(int *,uint64_t);
uint64_t *uint64Swap(uint64_t *,uint64_t);
float *floSwap(float *,uint64_t);
long int *lintSwap(long int *,uint64_t);
double *doSwap(double *,uint64_t);
double doOneSwap(double);
float floOneSwap(float);
uint32_t u32OneSwap(uint32_t);
int16_t *int16Swap(int16_t *,uint64_t);
double gaussian(double,double,double);
float randGauss(float,float);
double logNormal(double,double,double);
char *markChar(int,char *,char);
unsigned char *markUchar(int,unsigned char *,unsigned char);

unsigned char **uchChalloc(int,char *,int);
char **chChalloc(int,char *,int);
int **iIalloc(int,char *,int);
float **fFalloc(int,char *,int);
short int **shIalloc(int,char *,int);
double **dDalloc(int,char *,int);


void checkArguments(int,int,int,char *);
void TTIDY(void **,int); 

float singleMedian(float *,int);


