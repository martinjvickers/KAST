#include <math.h>

double d2(long long int ref[], long long int qry[], int nokmers);
double newd2(unordered_map<string, long long int> refmap, unordered_map<string, long long int> querymap);
double neweuler(unordered_map<string, long long int> refmap, unordered_map<string, long long int> querymap);
double euler(long long int ref[], long long int qry[], int nokmers);
double newd2s(unordered_map<string,thingy> qrymarkovthingy, unordered_map<string,thingy> refmarkovthingy);
double d2s(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[]);
double d2star(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[]);
