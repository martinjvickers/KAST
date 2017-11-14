#!/usr/bin/python
import os, sys, shutil, string
#import subprocess
import subprocess
import math
from operator import itemgetter, attrgetter
from optparse import OptionParser
from ALFSC_utils_v002 import *




#import matplotlib.pyplot as plt
def main( argv ) :

    parser = OptionParser()
    parser.add_option("-k", "--kmer", dest="klen", help="K-mer length [default: %default]", default=4 )
    parser.add_option("-p", "--prefix", dest="mylabel", help="Unique text for output prefix [default: %default]", default="kaog")
    parser.add_option("-q", "--query", dest="qfile", help="Query fasta file, multiple sequences [default: %default]", default="qry.fasta")
    parser.add_option("-r", "--reference", dest="rfile", help="Reference database fasta file, multiple sequences [default: %default]", default="ref.fasta")
    parser.add_option("-t", "--type", dest="type", help="Analysis type, one of: d2, d2s, d2star, d2s-opt, d2star-opt  [default: %default]", default="d2")
    parser.add_option("-c", "--revcomp", dest="revcomp", help="Include reverse complement of sequences [default: %default]", default='True')
    parser.add_option("-n", "--maxNumber", dest="maxcount", help="Maximum query sequences to analyse [default: %default]", default=1000)
    parser.add_option("-m", "--markov", dest="order", help="Markov order for d2s scoring type only [default: %default]", default=1)
    #parser.add_option("-f", "--topFeatures", dest="topBoxes", help="Only process boxes (features) with highest density [default: %default]", default=16384)
    #parser.add_option("-u", "--unitResolution", dest="resolution", help="Resolution, i.e. no of boxes in, or length of, unit side (set to power of 2) [default: %default]", default=9)

    (options, args) = parser.parse_args()
    if len(argv) <= 1 :
       parser.print_help()
       exit(0)

    print options
    print args

    if options.revcomp == 'True' : 
       doRC = 1
    else : 
       doRC = 0
    print "Setting doRC:",doRC
    kmerln = int(options.klen)
    myLabel = options.mylabel
    qfile = options.qfile
    rfile = options.rfile
    type = options.type
    maxCount = int(options.maxcount)
    markovOrder = int(options.order)

    #numOfHighestCounts = int(options.topBoxes)
    #nOfBoxesRoot = 2**int(options.klen)
    #unitLength = 2**(int(options.resolution))
    #numToKeep = numOfHighestCounts

    if int(options.klen) <= int(options.order) and type in ["d2s","d2star","d2s-opt","d2star-opt"]: 
       print "Exiting..."
       print "Markov order must be larger than kmer size otherwise d2S will not compute!\n",options.klen,options.order
       exit(0)

#    if int(options.klen) > int(options.resolution) :
#       print "Exiting"
#       print "Kmer length cannot be greater than unit side resolution (length)\n"
#       exit(0)

    output = myLabel+"."+str(kmerln)+"."+type
    if type in ["d2s","d2star","d2s-opt","d2star-opt"]:
       output = output+"."+str(markovOrder)
    print "Output prefix",output
    
    print "Reading sequences"
    seqs = {}
    maxSize = 1

    kmerSet = set([]) 
    seqs,QryKmerSet = processFasta2kmers( qfile, maxSize, kmerln, kmerSet )
    print "Read query fasta "+qfile+", number of seqs",len(seqs)
    kmerSet = set([])
    print "Number of kmers in QrykmerSet",len(QryKmerSet),"from",len(seqs),"sequences"

    refs,RefKmerSet = processFasta2kmers( rfile, maxSize, kmerln, kmerSet )
    print "Read refs fasta "+rfile+", number of seqs",len(refs)
    kmerSet = set([])
    print "Number of kmers in RefkmerSet",len(RefKmerSet),"from",len(refs),"sequences"

    unionKmerSet = RefKmerSet.union(QryKmerSet)
    print "Number of kmers in kmerSet",len(unionKmerSet),"max for kmer length",kmerln,"is",pow(4,kmerln),"difference",(pow(4,kmerln)-len(unionKmerSet))

#    seqs = processFasta( qfile, maxSize  )
#    print "Read query fasta "+qfile+", number of seqs",len(seqs)
#    refs = processFasta( rfile, maxSize  )
#    print "Read refs fasta "+rfile+", number of seqs",len(refs)
#    data = refsFromFasta( refs, doRC, type, markovOrder, kmerln )

    if type in ["d2s-opt","d2star-opt"]:
       data = refsFromFastaKS( refs, doRC, type, markovOrder, kmerln, unionKmerSet )
    else : 
       data = refsFromFasta( refs, doRC, type, markovOrder, kmerln )
    if "counts" in data : 
       allRefCounts = data["counts"]
       print "Got references: number "+str(len(allRefCounts))+" with kmer length",kmerln
    if "markov" in data : 
       print "Got references with markov order",
       allRefMarkov = data["markov"]  

    sizes = {}
    
    count = 1
    for id in seqs.iterkeys() :
       print "###",count,maxCount,"Doing QRY",id,"length",len(seqs[id])
       
       pnum = 1
       fragment = seqs[id]
       flen = len(fragment)
       size = flen
       sizes[id] = size
      
       print "--Fragment length: "+str(flen),id
       print "Calculating CGR matrix"       
       if doRC == 0 :
          allQseq = fragment
       else :
          rcfrag = reversecomplement(fragment)
          allQseq = fragment + "NNN" + rcfrag

       if type in ["kmer","d2","euler"] :
          scoresSorted = scoreCounts( allRefCounts, allQseq, kmerln, type ) 
       elif type in ["d2s","d2star"] : 
          print "Scoring with markov order:",markovOrder
          scoresSorted = scoreMarkov( allRefCounts, allRefMarkov, allQseq, kmerln, markovOrder, type )
       elif type in ["d2s-opt","d2star-opt"] :
          print "Scoring with markov order:",markovOrder
          scoresSorted = scoreMarkovOpt( allRefCounts, allRefMarkov, RefKmerSet, allQseq, kmerln, markovOrder, type )
       else :
           print "ERROR: Unknown scoring in main",type
           exit(1)
      
       print "Number of scores",len(scoresSorted) 
       print scoresSorted

       #continue
       outfile = output+".rank"
       f_o = open(outfile,"w")
       print "Outputing scores",outfile
       for i in range(len(scoresSorted)) :
          #print "RANK",i,scoresSorted[i],id,scoresSorted[i][1]
          CtgInfo = scoresSorted[i][1]
          #length = str(sizes[id])
          if CtgInfo.count("NODE") > 0 and CtgInfo.count("cov_") > 0 :
             vars = CtgInfo.split("_")
             coverage = vars[5]
             length = vars[3]
             winID = vars[7]
          elif CtgInfo.count("supercontig_supercontig") > 0 : 
             vars = CtgInfo.split(":")
             coverage = "n/a"
             length = "n/a"
          else :
             coverage = "n/a"
             length = "n/a"

          outStr = str(i)+"\t"+str(scoresSorted[i][0])+"\t"+length+"\t"+coverage+"\t"+winID+"\t"+CtgInfo+"\t"+id
          print "RANK",outStr
          f_o.write(outStr)
          f_o.write("\n")
       f_o.close()
       print "Written file",outfile

       count = count+1
       if count > maxCount : 
          print "BREAKING!",count
          break
   
    print "All done!"

main(sys.argv)

