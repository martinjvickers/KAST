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
#include <unordered_map>
#include <vector>
#include "common.h"

// Old Implementations
// Deprecated
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

double bray_curtis_distance(ModifyStringOptions options,
                            map<string, unsigned int> refcounts,
                            map<string, unsigned int> querycounts);

double normalised_google_distance(ModifyStringOptions options,
                                  map<string, unsigned int> refcounts,
                                  map<string, unsigned int> querycounts);

double chebyshev(ModifyStringOptions options,
                 map<string, unsigned int> refcounts,
                 map<string, unsigned int> querycounts);

double manhattan(ModifyStringOptions options,
                 map<string, unsigned int> refcounts,
                 map<string, unsigned int> querycounts);

double d2(ModifyStringOptions options, map<string, unsigned int> refcounts,
          map<string, unsigned int> querycounts);

double dai(ModifyStringOptions options, map<string, bool> ourkmers,
              map<string, unsigned int> refcounts,
              map<string, double> refmarkov,
              map<string, unsigned int> querycounts,
              map<string, double> querymarkov);


// New Template Implementations
template <typename TAlphabet>
double d2s(vector<String<TAlphabet>> ourkmers,
           map<String<TAlphabet>, unsigned int> refcounts, map<String<TAlphabet>, double> refmarkov,
           map<String<TAlphabet>, unsigned int> querycounts,
           map<String<TAlphabet>, double> querymarkov)
{
   double score = 0.0;
   double D2S = 0.0;
   double sum1 = 0.0;
   double sum2 = 0.0;

   unsigned int qN = 0;
   unsigned int rN = 0;

   for(auto p : querycounts)
   {
      qN = qN + p.second;
   }

   for(auto r : refcounts)
   {
      rN = rN + r.second;
   }

   for(auto p : ourkmers)
   {
      double qC = querycounts[p];
      double qP = querymarkov[p];
      double rC = refcounts[p];
      double rP = refmarkov[p];

      double qCt = qC - (qN*qP);
      double rCt = rC - (rN*rP);

      double dist = sqrt(qCt*qCt + rCt*rCt);
      if(dist == 0)
         dist = 1;
      D2S = D2S + (qCt*rCt / dist);
      sum1 = sum1 + (qCt*qCt / dist);
      sum2 = sum2 + (rCt*rCt / dist);
   }

   score = 0.5 * (1 - D2S/( sqrt(sum1)*sqrt(sum2) ));
   return score;
};

template <typename TAlphabet>
double d2star(vector<String<TAlphabet>> &ourkmers,
           map<String<TAlphabet>, unsigned int> &refcounts, map<String<TAlphabet>, double> &refmarkov,
           map<String<TAlphabet>, unsigned int> &querycounts,
           map<String<TAlphabet>, double> &querymarkov)
{
   double D_2Star = 0.0;
   double tempQ = 0.0;
   double tempR = 0.0;

   unsigned int qN = 0;
   unsigned int rN = 0;

   for(auto p : querycounts)
   {
      qN = qN + p.second;
   }

   for(auto r : refcounts)
   {
      rN = rN + r.second;
   }

   for(auto p : ourkmers)
   {
      double qC = querycounts[p];
      double qP = querymarkov[p];
      double rC = refcounts[p];
      double rP = refmarkov[p];

      double temp1 = qN * qP;
      double temp2 = rN * rP;
      double cQi_bar = qC - temp1;
      double cRi_bar = rC - temp2;
      double temp3 = sqrt(temp1*temp2);

      if(temp1 == 0)
      {
         temp1 = 1.0;
         temp3 = 1.0;
      }
      if(temp2 == 0)
      {
         temp2 = 1.0;
         temp3 = 1.0;
      }

      D_2Star += cQi_bar * cRi_bar / temp3;
      tempQ += cQi_bar * cQi_bar / temp1;
      tempR += cRi_bar * cRi_bar / temp2;
   }
   double temp = D_2Star / (sqrt(tempQ) * sqrt(tempR));
   return 0.5 * (1 - temp);
};

template <typename TAlphabet>
double hao(vector<String<TAlphabet>> ourkmers,
           map<String<TAlphabet>, unsigned int> refcounts, map<String<TAlphabet>, double> refmarkov,
           map<String<TAlphabet>, unsigned int> querycounts,
           map<String<TAlphabet>, double> querymarkov)
{
   int qN = 0;
   int rN = 0;

   for(auto p : querycounts)
   {
      qN = qN + p.second;
   }

   for(auto r : refcounts)
   {
      rN = rN + r.second;
   }

   double tempX = 0.0;
   double tempY = 0.0;
   double tempXY = 0.0;
   double d_Hao = 0.0;

   for(auto p : ourkmers)
   {
      double qC = querycounts[p];
      double qP = querymarkov[p];
      double rC = refcounts[p];
      double rP = refmarkov[p];

      double fQi = qC / (double)qN;
      double fRi = rC / (double)rN;
      double temp1, temp2;
      if(qP == 0)
         temp1 = -1;
      else
         temp1 = (fQi / qP) - 1.0;

      if(rP == 0)
         temp2 = -1;
      else
         temp2 = (fRi / rP) - 1.0;

      tempXY += temp1 * temp2;
      tempX += temp1 * temp1;
      tempY += temp2 * temp2;
   }
   double temp = tempXY / (sqrt(tempX) * sqrt(tempY));
   return (1 - temp) / 2;
};

template <typename TAlphabet>
double euler(map<String<TAlphabet>, unsigned int> &refcounts,
             map<String<TAlphabet>, unsigned int> &querycounts)
{
   double score = 0.0;
   long long unsigned int rN = 0;
   long long unsigned int qN = 0;

   for(pair<String<TAlphabet>, unsigned int> p: refcounts)
      rN = rN + p.second;
   for(pair<String<TAlphabet>, unsigned int> p: querycounts)
      qN = qN + p.second;

   double rN_new = (double)rN;
   double qN_new = (double)qN;

   // create a unified map
   map<String<TAlphabet>, unsigned int> ourkmers;
   ourkmers = refcounts;
   ourkmers.insert(querycounts.begin(),querycounts.end());

   // sum 
   for(pair<String<TAlphabet>, unsigned int> p: ourkmers)
   {
      double rF = refcounts[p.first] / (double)rN;
      double qF = querycounts[p.first] / (double)qN;
      score = score + (pow((rF - qF), 2));
   }

   return pow(score, 0.5);
};

template <typename TAlphabet>
double d2(map<String<TAlphabet>, unsigned int> refcounts,
          map<String<TAlphabet>, unsigned int> querycounts)
{
   long long unsigned int sumqCrC = 0;
   long long unsigned int sumqC2 = 0;
   long long unsigned int sumrC2 = 0;

   //create a unified map
   map<String<TAlphabet>, unsigned int> ourkmers;
   ourkmers = refcounts;
   ourkmers.insert(querycounts.begin(),querycounts.end());

   for(pair<String<TAlphabet>, unsigned int> p: ourkmers)
   {
      sumqCrC += ((long long unsigned int)refcounts[p.first] * (long long unsigned int)querycounts[p.first]);
      sumqC2 += ((long long unsigned int)querycounts[p.first] * (long long unsigned int)querycounts[p.first]);
      sumrC2 += ((long long unsigned int)refcounts[p.first] * (long long unsigned int)refcounts[p.first]);
    }

    long double score = sumqCrC / ((double)sqrt(sumqC2) * (double)sqrt(sumrC2));
    return 0.5*(1-score);
};

// New Template Implementation
template <typename TAlphabet>
double manhattan(map<String<TAlphabet>, unsigned int> refcounts,
                 map<String<TAlphabet>, unsigned int> querycounts)
{
   double score = 0.0;
   long long unsigned int rN = 0;
   long long unsigned int qN = 0;

   for(pair<String<TAlphabet>, unsigned int> p: refcounts)
      rN = rN + p.second;

   for(pair<String<TAlphabet>, unsigned int> p: querycounts)
      qN = qN + p.second;

   //create a unified map
   map<String<TAlphabet>, unsigned int> ourkmers;
   ourkmers = refcounts;
   ourkmers.insert(querycounts.begin(), querycounts.end());

   for(pair<String<TAlphabet>, unsigned int> p: ourkmers)
   {
      double rF = refcounts[p.first] / (double)rN;
      double qF = querycounts[p.first] / (double)qN;
      score = score + (abs(rF - qF));
   }
   return score;
};

template <typename TAlphabet>
double chebyshev(map<String<TAlphabet>, unsigned int> refcounts,
                 map<String<TAlphabet>, unsigned int> querycounts)
{
   double score = 0.0;
   long long unsigned int rN = 0;
   long long unsigned int qN = 0;
   double temp = 0.0;

   for(pair<String<TAlphabet>, unsigned int> p: refcounts)
      rN = rN + p.second;

   for(pair<String<TAlphabet>, unsigned int> p: querycounts)
      qN = qN + p.second;

   //create a unified map
   map<String<TAlphabet>, unsigned int> ourkmers;
   ourkmers = refcounts;
   ourkmers.insert(querycounts.begin(),querycounts.end());

   for(pair<String<TAlphabet>, unsigned int> p: ourkmers)
   {
      double rF = refcounts[p.first] / (double)rN;
      double qF = querycounts[p.first] / (double)qN;
      temp = abs(rF - qF);
      if(temp > score)
         score = temp;
   }
   return score;
};

template <typename TAlphabet>
double normalised_google_distance(map<String<TAlphabet>, unsigned int> refcounts,
                                  map<String<TAlphabet>, unsigned int> querycounts)
{
   double score = 0.0;
   double sumqC = 0.0;
   double sumrC = 0.0;

   map<String<TAlphabet>, unsigned int> min_qr;

   for(pair<String<TAlphabet>, unsigned int> p: querycounts)
      sumqC += p.second;

   for(pair<String<TAlphabet>, unsigned int> p: refcounts)
   {
      sumrC += p.second;
      if(querycounts.find(p.first) != querycounts.end())
      {
         if(querycounts[p.first] < p.second)
            min_qr[p.first] = querycounts[p.first];
         else
            min_qr[p.first] = p.second;
      }
   }

   double sum_min_qr = 0.0;
   for(pair<String<TAlphabet>, unsigned int> p: min_qr)
      sum_min_qr += min_qr[p.first];

   double sum_max, sum_min;

   if(sumqC > sumrC)
   {
      sum_max = sumqC;
      sum_min = sumrC;
   }
   else
   {
      sum_max = sumrC;
      sum_min = sumqC;
   }

   double sum_all = sumqC + sumrC;

   return (sum_max - sum_min_qr) / (sum_all - sum_min);
};

template <typename TAlphabet>
double bray_curtis_distance(map<String<TAlphabet>, unsigned int> refcounts,
                            map<String<TAlphabet>, unsigned int> querycounts)
{
   double sumMinus = 0.0;
   double sumPlus = 0.0;

   //create a unified map
   map<String<TAlphabet>, unsigned int> ourkmers;
   ourkmers = refcounts;
   ourkmers.insert(querycounts.begin(), querycounts.end());

   for(pair<String<TAlphabet>, unsigned int> p: ourkmers)
   {
      int qC = querycounts[p.first];
      int rC = refcounts[p.first];
      sumMinus = sumMinus + abs(qC - rC);
      sumPlus = sumPlus + abs(qC + rC);
   }

   return (double)sumMinus/(double)sumPlus;
};

template <typename TAlphabet>
double dai(vector<String<TAlphabet>> ourkmers,
           map<String<TAlphabet>, unsigned int> refcounts, map<String<TAlphabet>, double> refmarkov,
           map<String<TAlphabet>, unsigned int> querycounts,
           map<String<TAlphabet>, double> querymarkov)
{
   double score = 0.0;
   double S2kr = 0.0;
   int n = 0;
   double tempX = 0.0;
   double tempY = 0.0;

   int qN = 0;
   int rN = 0;

   for(auto p : querycounts)
   {
      qN = qN + p.second;
   }

   for(auto r : refcounts)
   {
      rN = rN + r.second;
   }

   for(auto p: ourkmers)
   {
      n++;
      double rC = refcounts[p];
      double qC = querycounts[p];
      double rP = refmarkov[p];
      double qP = querymarkov[p];

      double fRi = rC / (float)rN;
      double fQi = qC / (float)qN;
      double r_sigma = fRi * rP;
      double q_sigma = fQi * qP;
      double temp1 = 0.0;
      double temp2 = 0.0;

      // I added this if statement here. I'm not sure if this breaks
      // something, but basically we only do the calculation on kmers
      // that exist in both sequences. This d2tools implementions
      // simply bombs out if this situation is encountered
      if(rC != 0 && qC != 0)
      {
         if(r_sigma != q_sigma)
         {
            double temp = (float)(r_sigma+q_sigma);
            double temp3 = (2 * r_sigma) / (float)temp;
            double temp4 = (2 * q_sigma) / (float)temp;
            temp1 = r_sigma*log(temp3);
            temp2 = q_sigma*log(temp4);
         }
         tempX += temp1;
         tempY += temp2;
      }
   }

   if(tempX == 0.0 && tempY == 0.0)
      return 0.0;
   else
      return (tempX+tempY)/(n+2*log(2.0));

   return score;
};


#endif
