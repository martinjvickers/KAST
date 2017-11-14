#!/usr/bin/python
#import os, sys, shutil, string
import string
#import subprocess
import math
#import random
from decimal import *
from operator import itemgetter, attrgetter
import time

#import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.pyplot as plt


def reversecomplement(sequence):

    valid = ["A","C","G","T","N","c","g","t","a","n"]
    complement = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "a":"t", "t":"a", "c":"g", "g":"c", "n":"n"}

    sequence = sequence[::-1]

    reverse_complement_sequence = ""
    for letter in sequence:
        if letter not in valid :
           letter = "N"
        reverse_complement_sequence += complement[letter]

    return reverse_complement_sequence


def countKmersWindow (  allSeq, klen ) :

   p1 = 0
   p2 = p1+klen
   kmerCts = {}
   kmers = set([])

   if p2 >= len(allSeq) :
      print "Kmer length longer than sequence."
      print "No kmers!!"
 
   while p2 <= len(allSeq) : 

      kmer = allSeq[p1:p2]
      if "N" not in kmer :
         
         if kmer not in kmerCts : 
            kmerCts[kmer] = 0
            kmers.add(kmer)
         kmerCts[kmer] += 1

         p1 += 1
         p2 = p1+klen

      else :
     
        skip = kmer.find("N")
        p1 += skip+1
        p2 = p1+klen

   return kmerCts,kmers



def refsFromFastaKS( refs, doRC, type, markovOrder, klen, kmerSet ) :

   allCounts = {}
   allMarProb = {}
   wordToProb = {}

   #kmers = []
   #for kmer in kmerSet : 
   #   kmers.append(kmer)

   countAllKmers = 1
   num = 0
   for fid in refs.iterkeys() :

      print "REF Doing",num,fid,refs
      #id = id2prefix( fid )
      #fid = id
      print fid
      sequence = refs[fid]
      print "REF ID: "+fid+" length "+str(len(sequence))

      if doRC == 0 :
         allSeq = sequence
      else :
         revcomp = reversecomplement(sequence)
         allSeq = sequence+"NNN"+revcomp

      kcounts = countKmers( kmerSet, allSeq, klen )
      allCounts[fid] = kcounts
      #if countAllKmers == 1 :
      #i   refKmers = defineKmers( klen )
      #   kcounts = countKmers( refKmers, allSeq, klen )
      #   allCounts[id] = kcounts
      #else :
      #   kmerCts,rKmers = countKmersWindow ( allSeq, klen )
      #   allCounts[id] = kmerCts
         #print len(refKmerSet),pow(4,klen)
         #if len(refKmers) == pow(4,klen) :
         #   print "Got full complement of reference kmers"
         #else :
         #   for kmer in kmerCts.iterkeys() :
         #      if kmer not in refKmers :
         #         refKmers.append(kmer)
         #print len(refKmerSet),pow(4,klen)
         #for kmer in refKmerSet :
         #   if kmer in kmerCts :
         #      #print kmer,"n/a",kcounts[kmer]
         #   #else:
         #      print kmer,kmerCts[kmer],kcounts[kmer]
         #exit()

      if type.count("d2s") > 0 :
         allMarProb[fid] = getMarkovProbs( markovOrder, kmerSet, allSeq  )

      num = num+1

   print "Calculated counts!",num
   print "Number of counts:",len(allCounts)
   print "Number of markov probabilities",len(allMarProb),"for order",markovOrder

   data = {}
   data["counts"] = allCounts
   if type.count("d2s") > 0 :
      data["markov"] = allMarProb
   #data["refKmers"] = refKmers

   print "Number of sets in returned data",len(data)

   return data



def refsFromFasta( refs, doRC, type, markovOrder, klen ) :

   allCounts = {}
   allMarProb = {}
   wordToProb = {}
   refKmerSet = []

   countAllKmers = 1
   num = 0
   for fid in refs.iterkeys() :

      print "REF Doing",num,fid
      #id = id2prefix( fid )
      #id = fid
      sequence = refs[fid]
      print "REF ID: "+fid+" length "+str(len(sequence))

      if doRC == 0 :
         allSeq = sequence
      else : 
         revcomp = reversecomplement(sequence)
         allSeq = sequence+"NNN"+revcomp

      #if countAllKmers == 1 : 
      refKmers = defineKmers( klen )
      kcounts = countKmers( refKmers, allSeq, klen )
      allCounts[fid] = kcounts
      #else :
      #   kmerCts = countKmersWindow ( allSeq, klen ) 
      #   allCounts[id] = kmerCts
      #   #print len(refKmerSet),pow(4,klen)
      #   if len(refKmers) == pow(4,klen) :
      #      print "Got full complement of reference kmers"
      #   else : 
      ##      for kmer in kmerCts.iterkeys() : 
       #        if kmer not in refKmers : 
       #           refKmers.append(kmer)     
         #print len(refKmerSet),pow(4,klen)
         #for kmer in refKmerSet :
         #   if kmer in kmerCts :
         #      #print kmer,"n/a",kcounts[kmer]
         #   #else:
         #      print kmer,kmerCts[kmer],kcounts[kmer]
         #exit()

      if type in ["d2s","d2star"] : 
         allMarProb[fid] = getMarkovProbs( markovOrder, refKmers, allSeq  )

      num = num+1

   print "Calculated counts!",num
   print "Number of counts:",len(allCounts)
   print "Number of markov probabilities",len(allMarProb),"for order",markovOrder

   data = {}
   data["counts"] = allCounts
   if type in ["d2s","d2star"] :
      data["markov"] = allMarProb
   data["refKmers"] = refKmers

   print "Number of sets in returned data",len(data)   

   return data

def getMarkov( markovOrder, kmers, sequence ) :

   klen = markovOrder

   # get probability of markov order kmers
   mpKmers = defineKmers( klen )
   print mpKmers
   mpKcount = countKmers( mpKmers, sequence, klen )
   print mpKcount

   num = 0.0
   for k in mpKcount.iterkeys() :
      num += mpKcount[k]

   print "Probabilities for markov kmers, number:",num,"Markov order",klen
   mpKfreq = {}
   tot = 0.0
   for k in mpKcount.iterkeys() :
      mpKfreq[k] = float(mpKcount[k]) / num
      #print k, mpKcount[k], mpKfreq[k]
      tot += mpKfreq[k]
   print "Probability sum:",tot

   return mpKfreq


def getMarkovProbs( markovOrder, kmers, sequence  ) :

   mpKfreq = getMarkov( markovOrder, kmers, sequence )
   #print mpKfreq
   # now get probability of kmers using the mpKmers
   print "Probabilities for kmers from markov model"
   kprob = {}
   tot = 0.0
   for kmer in kmers :
      prob = 1.0
      for i in range(len(kmer)) :
         j = i+markovOrder
         word = kmer[i:j]
         prob = prob * mpKfreq[word]
         #print "mN",word,kmer,prob,mpKfreq[word],markovOrder
         if j == len(kmer) :
            break
      kprob[kmer] = prob
      #print kmer,prob
      tot += prob
   print "Probability sum:",tot
 
   return kprob  

#def getMarkovProbs( markovOrder, kmers, sequence  ) :
#
#   if markovOrder >= 0 :
#      mpKfreq0 = getMarkovFreq( 0, kmers, sequence )
#   if markovOrder >= 1 :
#      mpKfreq1 = getMarkovFreq( 1, kmers, sequence )
#   if markovOrder >= 2 : 
#      print "WARNING: Ignoring Markov Orders of 2 and above!"
#
#   # now get probability of kmers using the mpKmers
#   print "Probabilities for kmers from markov model"
#   kprob = {}
#   tot = 0.0
#   for kmer in kmers :
#      prob = 1.0
#      for i in range(len(kmer)) :
#         if i == 0 and markovOrder > 0: 
#            bp = kmer[-1]
#            prob *= mpKfreq0[bp]
#            print "m0",bp,kmer,prob,mpKfreq0[bp]
##         j = i+markovOrder+1
#         word = kmer[i:j]
#         if markovOrder == 0 : 
#            prob = prob * mpKfreq0[word]
#            print "mN",word,kmer,prob,mpKfreq0[word],markovOrder
#         elif markovOrder == 1 :
#            prob = prob * mpKfreq1[word]
#            print "mN",word,kmer,prob,mpKfreq1[word],markovOrder
##         if j == len(kmer) :
#            break
#      kprob[kmer] = prob
#      print kmer,prob
#      tot += prob
#   print "Probability sum:",tot
#
#   return kprob
#

def countKmers( kmers, sequence, klen ) :

   kmerCount = {}

   for kmer in kmers :
      kmerCount[kmer] = 0

   sequence = sequence.upper()
   seqL = []
   for char in list(sequence) : 
      if char not in ["A","C","G","T","N"] : 
         char = "N"
      seqL.append(char)

   max = len(seqL) - klen
   i = 0
   while i <= max :
      j = i + klen
      wordL = seqL[i:j]
      if "N" in wordL : 
         i = i+1
         continue
      word = "".join(wordL)
      kmerCount[word] = kmerCount[word] + 1
      i = i+1

   #for word in kmerCount.iterkeys() : 
   #   print "COUNT",word,kmerCount[word]

   #print "countKmers: num contigs in scaffold:",num
   return kmerCount

def scoreMarkovOpt( allRefCounts, allRefMarkov, refKmerSet, allQseq, klen, markovOrder, type ) :

   qcounts,qKmers = countKmersWindow ( allQseq, klen )
   kmerSet = refKmerSet.union(qKmers)
   qmarkovProb = getMarkovProbs( markovOrder, kmerSet, allQseq )
   print "kmerSet:",len(kmerSet),"versus",len(allQseq)
   print allQseq
   scores = []
   for ref in allRefCounts.iterkeys() :

      refcounts = allRefCounts[ref]
      rmarkovProb = allRefMarkov[ref]

      if type == "d2s-opt" or type == "d2star-opt" :

         D2Star = D2S = sum1 = sum2 = score = 0.0
         qN = rN = 0.0

         for kmer in kmerSet : 
            if kmer not in qcounts :
               qcounts[kmer] = 0

            if kmer not in refcounts :
               refcounts[kmer] = 0

            qN = qN + qcounts[kmer]
            rN = rN + refcounts[kmer]
            print kmer,"Q",qcounts[kmer],qN,"R",refcounts[kmer],rN

         for kmer in kmerSet :

            if qcounts[kmer] == 0 and refcounts[kmer] == 0 : 
               continue

            qC = float(qcounts[kmer])
            rC = float(refcounts[kmer])
            qP = qmarkovProb[kmer]
            rP = rmarkovProb[kmer]

            q_np = qN*qP
            r_np = rN*rP
            qCt = qC - q_np
            rCt = rC - r_np
            #print kmer,"Q",qCt,"from",qC,":",qN,qP,(qN*qP)
            #print kmer,"R",rCt,"from",rC,":",rN,rP,(rN*rP)

            if type == "d2s-opt" :
               dist = math.sqrt(qCt*qCt + rCt*rCt)
               if dist == 0.0 : 
                  print "WARNING: division by 0, continuing"
                  continue
               D2S = D2S + (qCt*rCt / dist)
               sum1 = sum1 + (qCt*qCt / dist)
               sum2 = sum2 + (rCt*rCt / dist)

            if type == "d2star-opt" :
               dist = math.sqrt(q_np*r_np)
               if dist == 0.0 :
                  print "WARNING: division by 0, continuing"
                  continue
               D2Star = D2Star + ((qCt*rCt) / dist)
               sum1 = sum1 + (qCt*qCt / q_np )
               sum2 = sum2 + (rCt*rCt / r_np )

         if type == "d2s-opt" :   
            #print "D2S",D2S
            score = 0.5 * (1 - D2S/( math.sqrt(sum1)*math.sqrt(sum2) ))

         if type == "d2star-opt" : 
            #print "D2Star",D2Star
            score = 0.5 * (1 - D2Star/( math.sqrt(sum1)*math.sqrt(sum2) ))

      else :
         print "ERROR: unknown scoring type in markov scoring:",type
         exit(1)

      if score == 999 :
         print "WARNING: could not score:",ref
         continue

      scores.append([score,ref])

   print "Sorting scores for all references, no. of scores:",len(scores)
   scoresSorted = sorted(scores, key=itemgetter(0))
   return scoresSorted


def scoreMarkov( allRefCounts, allRefMarkov, allQseq, klen, markovOrder, type ) :
   qKmers = defineKmers( klen )
   qcounts = countKmers( qKmers, allQseq, klen )
   qProbs = getMarkovProbs( markovOrder, qKmers, allQseq )
   #print "HERERERER"
   #print qcounts
   #print qProbs
   #exit()

   scores = []

   counter = 0

   for ref in allRefCounts.iterkeys() :

      refcounts = allRefCounts[ref]
      refprobs = allRefMarkov[ref]

      if type == "d2s" :
	 start = time.clock()
         score = d2s( qcounts, refcounts, qProbs, refprobs ) 
	 print "d2 Time taken", time.clock() - start
      elif type == "d2star" :
         score = d2star( qcounts, refcounts, qProbs, refprobs )
      else :
         print "ERROR: unknown scoring type in markov scoring:",type
         exit(1)

      if score == 999 : 
         print "WARNING: could not score:",ref
         continue

      #ref = "dummy/"+ref
      scores.append([score,ref])
      counter = counter + 1

   print "Sorting scores for all references, no. of scores:",len(scores)
   scoresSorted = sorted(scores, key=itemgetter(0))
   #scoresSorted.reverse()
   return scoresSorted


def scoreCounts( allRefCounts, allQseq, klen, type ) :

   qKmers = defineKmers( klen )
   qcounts = countKmers( qKmers, allQseq, klen )

   start = time.clock()

   scores = []
   for ref in allRefCounts.iterkeys() :

      refcounts = allRefCounts[ref]
      if type == "d2" :
         score = d2( qcounts, refcounts )
      elif type in ["euler","kmer"] :
         score = euler( qcounts, refcounts )
      else :
         print "ERROR: unknown scoring type:",type
         exit(1)

      #ref = "dummy/"+ref
      scores.append([score,ref])

   
   print "Sorting scores for all references, no. of scores:",len(scores)
   print "Time taken", time.clock() - start
   scoresSorted = sorted(scores, key=itemgetter(0))
   return scoresSorted


def d2star(  qcounts, refcounts, qmarkovProb, rmarkovProb ) :

   print "Doing D2Star"
   score = 0.0
   D2Star= 0.0
   sum1 = 0.0
   sum2 = 0.0

   rN = 0.0
   qN = 0.0
   for kmer in qcounts.iterkeys() :
      qC = qcounts[kmer]
      rC = refcounts[kmer]
      qN = qN + qC
      rN = rN + rC

   for kmer in qcounts.iterkeys() :
      qC = qcounts[kmer]
      rC = refcounts[kmer]
      qP = qmarkovProb[kmer]
      rP = rmarkovProb[kmer]

      q_np = qN*qP
      r_np = rN*rP

      qCt = qC - q_np
      rCt = rC - r_np
      #print kmer,"qCt","rCt",qCt,qC,qP,rCt,rC,rP

      #BUG???
      #dist = math.sqrt(q_np) * math.sqrt(r_np)
      dist = math.sqrt(q_np*r_np)
      #Should be: dist = math.sqrt(q_np*r_np)
      if dist == 0.0 :
         print "WARNING: Zero division!!"
         print "qC,qN,qP and qCt:",qC,qN,qP,qCt
         print "rC,rN,rP and rCt:",rC,rN,rP,rCt
         print "Adding 0.0 to sums / continuing"
         continue

      D2Star= D2Star + ((qCt*rCt) / dist)
      sum1 = sum1 + (qCt*qCt / q_np )
      sum2 = sum2 + (rCt*rCt / r_np )

   score = 0.5 * (1 - D2Star/( math.sqrt(sum1)*math.sqrt(sum2) ))
   #print "D2Star",D2Star,score
   return score

def d2s(  qcounts, refcounts, qmarkovProb, rmarkovProb ) :

   score = 0.0
   D2S = 0.0
   sum1 = 0.0
   sum2 = 0.0

   rN = 0.0
   qN = 0.0
   rtot = 0.0
   qtot = 0.0
   for kmer in qcounts.iterkeys() :
      qN = qN + qcounts[kmer]
      qtot = qtot + qmarkovProb[kmer]
      if kmer not in refcounts :
         print "Stopping, kmers not present in both ref and query. Try a smaller kmer?",kmer,len(refcounts),len(qcounts)
         exit()
      rN = rN + refcounts[kmer]
      rtot = rtot + rmarkovProb[kmer]

   #print "rN,qN",rN,qN
   #print "rtot,qtot",rtot,qtot

   #print qmarkovProb
   #print rmarkovProb
   counter = 0
   for kmer in qcounts.iterkeys() :

      qC = float(qcounts[kmer])
      rC = float(refcounts[kmer])
      qP = qmarkovProb[kmer]
      rP = rmarkovProb[kmer]

      qCt = qC - (qN*qP)
      rCt = rC - (rN*rP)
      #qCt = qC
      #rCt = rC
      #print kmer,"Q",qCt,"from",qC,":",qN,qP,(qN*qP)
      #print kmer,"R",rCt,"from",rC,":",rN,rP,(rN*rP)
      dist = math.sqrt(qCt*qCt + rCt*rCt)

      if dist == 0.0 :
         print "WARNING: Zero division!!"
         print "qC,qN,qP and qCt:",qC,qN,qP,qCt
         print "rC,rN,rP and rCt:",rC,rN,rP,rCt
         print "Adding 0.0 to sums / continuing"
         continue

      D2S = D2S + (qCt*rCt / dist)
      sum1 = sum1 + (qCt*qCt / dist)
      sum2 = sum2 + (rCt*rCt / dist)
      counter = counter + 1

   #print "D2S",D2S
   score = 0.5 * (1 - D2S/( math.sqrt(sum1)*math.sqrt(sum2) ))
   #print "D2S",D2S,score
   print "D2 loop", counter
   return score



def d2( qcounts, refcounts ) :

   sumqCrC = 0.0
   sumqC2 = 0.0
   sumrC2 = 0.0

   for kmer in qcounts.iterkeys() : 
      qC = qcounts[kmer]
      rC = refcounts[kmer]
      sumqCrC = sumqCrC + qC*rC
      sumqC2 = sumqC2 + qC*qC
      sumrC2 = sumrC2 + rC*rC

   score = sumqCrC / ( math.sqrt(sumqC2) * math.sqrt(sumrC2) )
   #print score
   d2 = 0.5*(1-score)

   return d2


def euler( qcounts, refcounts ) :

   score = 0.0
   rN = 0.0
   qN = 0.0

   for kmer in qcounts.iterkeys() : 
      qC = qcounts[kmer]
      rC = refcounts[kmer]
      qN = qN + qC
      rN = rN + rC
 
   for kmer in qcounts.iterkeys() :
      qC = qcounts[kmer]
      rC = refcounts[kmer]
      qF = qC / qN
      rF = rC / rN
      dist = (rF - qF)**2
      score = score + dist

   score = score**0.5

   return score


def distSq( x1,y1,x2,y2 ) :
   distSq = (x1-x2)**2 + (y1-y2)**2
   return distSq

def dist( x1,y1,x2,y2 ) :
   distSq = (x1-x2)**2 + (y1-y2)**2
   dist = math.sqrt(distSq)
   return dist


def getTaxInfo( gi, taxon ) :

   if gi not in taxon :
      taxinfo = "missing for "+gi
      tinfo = {}
   else :
      tinfo = taxon[gi]
      taxinfo = tinfo["spc"]+"::"+tinfo["gen"]+"::"+tinfo["fam"]
      taxinfo = taxinfo+"::"+tinfo["ord"]+"::"+tinfo["cls"]+"::"+tinfo["phy"]

   return tinfo,taxinfo


def sampledToTax( id, taxon ) :

   text = id.split("/")[1]
   text = text.split("_::")[0]
   gi = giFilenameToGI( text )
   taxRecord,taxStr = getTaxInfo( gi, taxon )
   tax = taxRecord['spc']+" || "+taxRecord['gen']+" || "+taxRecord['ord']+" || "+taxRecord['cls']

   return tax


def id2prefix( id ) :

   vars = id.split(",")[0].split()
   desc = vars[1:]
   gi = vars[0]
   gi = gi.replace("|","_")
   prefix = "_".join(desc)
   prefix = prefix.replace("/","-")
   prefix = gi+"_"+prefix

   return prefix


def getGItoNames( nfile ) :

   f_i = open(nfile,"r")
   name2gi = {}
   for line in f_i :
      vars = line.strip().split("\t")
      name2gi[vars[1].strip()] = vars[0].strip()
   f_i.close()
   return name2gi


def parseFileID( fileID ) :

   [p1,p2] = fileID.split("__")

   gi = giFilenameToGI ( p1 )

   name = p2.split(".")[0]

   return gi,name

def giFilenameToGI ( part ) :

   if part.count("gi_") > 0 :
      gi = part.replace("_","|")
      gi = gi.replace("NC|","NC_")
      gi = gi+"|"
      if gi.count("NZ|") > 0 :
         gi = gi.replace("NZ|","NZ_")

   return gi


def readTaxononmy( tfile ) :

   f_i = open(tfile,"r")

   taxon = {}
   for line in f_i :

     vars = line.strip().split(";;")

     gi = vars[0].strip()
     name = vars[1].strip()
   
     info = {}
     info["spc"] = "na"
     info["gen"] = "na"
     info["fam"] = "na"
     info["ord"] = "na"
     info["cls"] = "na"
     info["phy"] = "na"
     info["kng"] = "na"

     for item in vars :
  
        #print item 
        if item.count("species:") > 0 :
           info["spc"] = item.split(":")[1]#.split()[0]
        if item.count("genus:") > 0 :
           info["gen"] = item.split(":")[1]#.split()[0]
        if item.count("family:") > 0 :
           info["fam"] = item.split(":")[1]#.split()[0]
        if item.count("order:") > 0 :
           info["ord"] = item.split(":")[1]#.split()[0]
        if item.count("class:") > 0 :
           info["cls"] = item.split(":")[1]#.split()[0]
        if item.count("phylum:") > 0 :
           info["phy"] = item.split(":")[1].split()[0]
        if item.count("superkingdom:") > 0 :
           info["kng"] = item.split(":")[1]

     taxon[gi] = info      
     #print gi, name, taxon[gi]
    
   f_i.close()

   return taxon 


def processFasta( file, size ):

    f_in = open( file )

    Seqs = {}
    sequence = ""
    header = ""
    for line in f_in :

        line = line.strip()

        if line.count(">") == 1:

           if header != "":
              if len(sequence) >= size:
                 id = header.split(">")[1]
                 Seqs[id] = sequence.upper()

           header = line
           sequence = ""

        else:
           sequence = sequence+line

    if len(sequence) >= size :
       id = header.split(">")[1]
       Seqs[id] = sequence.upper()

    f_in.close()

    return Seqs


def processFasta2kmers( file, size, klen, kmerSet ):

    f_in = open( file )

    maxKmers = pow(4,klen)
    Seqs = {}
    sequence = ""
    header = ""
    for line in f_in :

        line = line.strip()

        if line.count(">") == 1:

           if header != "":
              if len(sequence) >= size:
                 id = header.split(">")[1]
                 Seqs[id] = sequence.upper()
                 if len(kmerSet) < maxKmers :
                    kmerSet = kmersWindow( Seqs[id], klen, kmerSet )

           header = line
           sequence = ""

        else:
           sequence = sequence+line

    if len(sequence) >= size :
       id = header.split(">")[1]
       Seqs[id] = sequence.upper()
       if len(kmerSet) < maxKmers  :
          kmerSet = kmersWindow( Seqs[id], klen, kmerSet )

    f_in.close()

    return Seqs,kmerSet

# used in the initial calculation of the universal kmer set
def kmersWindow (  seq, klen, kset ) :

   p1 = 0
   p2 = p1+klen
   maxKmers = pow(4,klen)

   if p2 >= len(seq) :
      print "Kmer length longer than sequence."
      print "No kmers!!"

   while p2 < len(seq) and len(kset) < maxKmers :

      kmer = seq[p1:p2]
      if "N" not in kmer :

         if kmer not in kset :
            rckmer = reversecomplement(kmer)
            kset.add(kmer)
            kset.add(rckmer)

         p1 += 1
         p2 = p1+klen

      else :

        skip = kmer.find("N")
        p1 += skip+1
        p2 = p1+klen

   return kset


def defineKmers( klen ) :

   base = ["A","C","G","T"]
   kmer = ""
   kmers = list(base)

   current = []
   current.append(kmer)

   for i in range(klen-1) :
      new = []
      for j in range(4) :
         for kmer in kmers :
            kmer = kmer+base[j]
            new.append(kmer)
      kmers = list(new)

   return kmers



