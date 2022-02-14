/*
KAST - Kmer Alignment-free Search Tool
Version 0.0.34
Written by Dr. Martin Vickers (martin.vickers@jic.ac.uk)

*/

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <math.h>
#include <string>
#include <thread>
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
   if(parseMask(options, options.effectiveLength) == 1)
      return 1;

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
   // Running in pairwise all mode
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
   // Running in interleaved mode
   else if(options.interleavedFileName != NULL)
   {
      if(options.sequenceType == "aa")
      {  
         interleaved(options, AminoAcid());
      }
      else if(options.sequenceType == "raa")
      {  
         interleaved(options, ReducedAminoAcidMurphy10());
      }
      else if(options.sequenceType == "dna")
      {  
         interleaved(options, Dna5());
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
