/*
KAST - Kmer Alignment-free Search Tool
Version 0.0.18
Written by Dr. Martin Vickers (martin.vickers@jic.ac.uk)

MIT License

Copyright (c) 2018 Martin James Vickers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>
//#include <seqan/store.h>
#include <string>
#include <thread>
//#include <map>
//#include <unordered_map>
//#include <vector>
//#include <seqan/alignment_free.h>
#include <seqan/reduced_aminoacid.h>
#include "common.h"
#include "utils.h"
#include "pairwise.h"
#include "search.h"

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv)
{
   // parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // parse the mask so we know the kmer size
   options.effectiveLength = options.klen;
   //if(parseMask(options, options.effectiveLength) == 1)
   //   return 1;

   // Running in pairwise mode
   if(options.pairwiseFileName != NULL && options.type != "all" && options.type != "new")
   {
      if(options.sequenceType == "aa")
      {
         pairwise_matrix(options, AminoAcid());
      }
      else if(options.sequenceType == "raa")
      {
         pairwise_matrix(options, ReducedAminoAcidMurphy10());
      }
      else if(options.sequenceType == "dna")
      {
         pairwise_matrix(options, Dna5());
      }
      else
      {
         // there is no other mode
         cerr << "Error: mode not found - " << options.sequenceType << endl;
         return 1;
      }
   }
   // Running in pairwise mode
   else if(options.pairwiseFileName != NULL && options.type == "all")
   {
      if(options.sequenceType == "aa")
      {
         pairwise_all_matrix(options, AminoAcid());
      }
      else if(options.sequenceType == "raa")
      {
         pairwise_all_matrix(options, ReducedAminoAcidMurphy10());
      }
      else if(options.sequenceType == "dna")
      {
         pairwise_all_matrix(options, Dna5());
      }
      else
      {
         // there is no other mode
         cerr << "Error: mode not found - " << options.sequenceType << endl;
         return 1;
      }
   }
   else if (options.referenceFileName != NULL && options.queryFileName != NULL)
   {
      if(options.sequenceType == "aa")
      {
         query_ref_search(options, AminoAcid());
      }
      else if(options.sequenceType == "raa")
      {
         query_ref_search(options, ReducedAminoAcidMurphy10());
      }
      else if(options.sequenceType == "dna")
      {
         query_ref_search(options, Dna5());
      }
      else
      {
         // there is no other mode
         cerr << "Error: mode not found - " << options.sequenceType << endl;
         return 1;
      }
   }

   return 0;
}
