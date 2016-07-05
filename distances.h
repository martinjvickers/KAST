#include <math.h>

double d2(long long int ref[], long long int qry[], int nokmers);
double euler(long long int ref[], long long int qry[], int nokmers);
double d2s(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[]);
double d2star(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[]);
