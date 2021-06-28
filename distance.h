#ifndef DIST_H
#define DIST_H

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <string>
#include "common.h"

double d2(String<unsigned> const & kmerCounts1,
          String<unsigned> const & kmerCounts2);

double euler(String<unsigned> const & kmerCounts1,
             String<unsigned> const & kmerCounts2);

double bray_curtis_distance(String<unsigned> const & kmerCounts1,
                            String<unsigned> const & kmerCounts2);

double normalised_google_distance(String<unsigned> const & kmerCounts1,
                                  String<unsigned> const & kmerCounts2);

double chebyshev(String<unsigned> const & kmerCounts1,
                 String<unsigned> const & kmerCounts2);

double canberra(String<unsigned> const & kmerCounts1,
                String<unsigned> const & kmerCounts2);

double normalised_canberra(String<unsigned> const & kmerCounts1,
                           String<unsigned> const & kmerCounts2);

double manhattan(String<unsigned> const & kmerCounts1,
                 String<unsigned> const & kmerCounts2);

double d2s(String<unsigned> const & kmerCounts1,
           String<unsigned> const & kmerCounts2,
           String<double> const & markovCounts1,
           String<double> const & markovCounts2);

double d2star(String<unsigned> const & kmerCounts1,
              String<unsigned> const & kmerCounts2,
              String<double> const & markovCounts1,
              String<double> const & markovCounts2);

double dai(String<unsigned> const & kmerCounts1,
           String<unsigned> const & kmerCounts2,
           String<double> const & markovCounts1,
           String<double> const & markovCounts2);

double hao(String<unsigned> const & kmerCounts1,
           String<unsigned> const & kmerCounts2,
           String<double> const & markovCounts1,
           String<double> const & markovCounts2);

#endif
