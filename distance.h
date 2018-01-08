#ifndef DIST_H
#define DIST_H

#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <string>
#include <thread>
#include <map>
#include <vector>
//#include <gmp.h>
#include "common.h"

double d2s(ModifyStringOptions options, map<string, bool> ourkmers, 
           map<string, unsigned int> refcounts, map<string, double> refmarkov, 
           map<string, unsigned int> querycounts, 
           map<string, double> querymarkov);
double d2star(ModifyStringOptions options, map<string, bool> ourkmers, 
              map<string, unsigned int> refcounts, 
              map<string, double> refmarkov, 
              map<string, unsigned int> querycounts, 
              map<string, double> querymarkov);
double hao(ModifyStringOptions options, map<string, bool> ourkmers, 
           map<string, unsigned int> refcounts, map<string, double> refmarkov, 
           map<string, unsigned int> querycounts, 
           map<string, double> querymarkov);
double euler(ModifyStringOptions options, map<string, unsigned int> refcounts, 
             map<string, unsigned int> querycounts);
double d2(ModifyStringOptions options, map<string, unsigned int> refcounts, 
          map<string, unsigned int> querycounts);
double manhattan(ModifyStringOptions options, 
                 map<string, unsigned int> refcounts, 
                 map<string, unsigned int> querycounts);
double chebyshev(ModifyStringOptions options, 
                 map<string, unsigned int> refcounts, 
                 map<string, unsigned int> querycounts);
double normalised_google_distance(ModifyStringOptions options, 
                                  map<string, unsigned int> refcounts, 
                                  map<string, unsigned int> querycounts);
double bray_curtis_distance(ModifyStringOptions options, 
                            map<string, unsigned int> refcounts, 
                            map<string, unsigned int> querycounts);
double dai(ModifyStringOptions options, map<string, bool> ourkmers, 
              map<string, unsigned int> refcounts, 
              map<string, double> refmarkov, 
              map<string, unsigned int> querycounts, 
              map<string, double> querymarkov);
#endif
