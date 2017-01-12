#include <unordered_map>

//calculate the counts
//BUG!!! I need a check that if the klen is longer than the sequence, it should fall out gracefully
count_obj count(int klen, IupacString sequence)
{
        int total = 0;
        unordered_map<string, long long int> map;

        /*
        if(klen==1)
        {
                cout << "testing count k=1 " << endl;
        }
        */

        //iterate over the sequence
        for(int i = 0; i <= length(sequence)-klen; i++)
        {
                //get our kmer
                string kmer;
                assign(kmer,infix(sequence, i, i+klen));

                //need to drop if there is an N in it
                size_t found = kmer.find("N");

                if(found > kmer.size()){
                        long long int count = map[kmer];
                        map[kmer] = count + 1;
                /*      if(klen==1)
                        {
                                cout << kmer << " " << kmer.size() << " " << count << endl;
                        }
                */
                        total++;
                }
        }

        count_obj obj;
        obj.kmer_counts = map;
        obj.total = total;
        return obj;
}


double euler(Seq ref, Seq query, ModifyStringOptions options)
{
	double score = 0.0;
	int rN = ref.getTotalCounts(options.klen);
	int qN = query.getTotalCounts(options.klen);	

	//extract our maps
	unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
	unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

	//create a unified map
	unordered_map<string, long long int> ourkmers;
        ourkmers = refkmers;
	ourkmers.insert(querykmers.begin(),querykmers.end());

	for(pair<string, long long int> p: ourkmers)
        {
                double rF = refkmers[p.first] / (double)rN;
                double qF = querykmers[p.first] / (double)qN;
                score = score + (pow((rF - qF), 2));
        }

	return pow(score, 0.5);
}

double d2(Seq ref, Seq query, ModifyStringOptions options)
{
	double sumqCrC = 0.0;
	double sumqC2 = 0.0;
	double sumrC2 = 0.0;

	//extract our maps
        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

        //create a unified map
        unordered_map<string, long long int> ourkmers;
        ourkmers = refkmers;
        ourkmers.insert(querykmers.begin(),querykmers.end());

	for(pair<string, long long int> p: ourkmers)
        {
		sumqCrC = sumqCrC + (refkmers[p.first] * querykmers[p.first]);
                sumqC2 = sumqC2 + (querykmers[p.first] * querykmers[p.first]);
                sumrC2 = sumrC2 + (refkmers[p.first] * refkmers[p.first]);
        }
	
	double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}

double d2s(Seq ref, Seq query, ModifyStringOptions options)
{
	double score = 0.0;
        double D2Star = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;

	unordered_map<string,markov_dat> qrymarkov = query.getMarkov(options.klen, options.markovOrder);
	unordered_map<string,markov_dat> refmarkov = ref.getMarkov(options.klen, options.markovOrder);

	//create union of the two maps
        unordered_map<string, markov_dat> ourkmers;
        ourkmers = refmarkov;
        ourkmers.insert(qrymarkov.begin(),qrymarkov.end());

	for(pair<string, markov_dat> p: ourkmers)
        {
		double qC = qrymarkov[p.first].count;
                double rC = refmarkov[p.first].count;
                double qP = qrymarkov[p.first].prob;
                double rP = refmarkov[p.first].prob;

		double q_np = query.getTotalMarkov(options.klen, options.markovOrder) * qP;
		double r_np = ref.getTotalMarkov(options.klen, options.markovOrder) * rP;
		double qCt = (double)qC - q_np;
		double rCt = (double)rC - r_np;

		double dist = sqrt(( qCt * qCt )+( rCt * rCt ));

		D2Star = D2Star + ( ( qCt * rCt ) / dist );
		sum1 = sum1 + ( ( qCt * qCt ) / dist );
                sum2 = sum2 + ( ( rCt * rCt ) / dist );

		cout << "QRY: " << p.first << " " << qC << " " << qP << " " << q_np << " " << qCt << " " << dist << " " << D2Star << " " << sum1 << " " << sum2 << endl;
		cout << "REF: " << p.first << " " << rC << " " << rP << " " << r_np << " " << rCt << " " << dist << " " << D2Star << " " << sum1 << " " << sum2 << endl;
	}

	score = 0.5 * (1 - ( (D2Star) / (sqrt(sum1)*sqrt(sum2)) ) );
	return score;
}

double d2snew(Seq ref, Seq query, ModifyStringOptions options)
{
	int rN = ref.getTotalCounts(options.klen);
	int qN = query.getTotalCounts(options.klen);
	//double rN = ref.getSumProb(options.klen, options.markovOrder);
	//double qN = query.getSumProb(options.klen, options.markovOrder);
	double score = 0.0;
	double D2S = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	
	//get my markov maps
	unordered_map<string,markov_dat> qrymarkov = query.getMarkov(options.klen, options.markovOrder);
	unordered_map<string,markov_dat> refmarkov = ref.getMarkov(options.klen, options.markovOrder);
	
        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

	//create union of the two maps
        unordered_map<string, markov_dat> ourkmers;
        ourkmers = refmarkov;
        ourkmers.insert(qrymarkov.begin(),qrymarkov.end());

	for(pair<string, markov_dat> p: ourkmers)
        {

//		double qC = qrymarkov[p.first].count;
//                double rC = refmarkov[p.first].count;
		double qC = querykmers[p.first];
		double qP;
		if(qC==0.0)
		{
			qP = query.missing(p.first,options.markovOrder);
		} else {
			qP = qrymarkov[p.first].prob;
		}
		
		double rC = refkmers[p.first];
		double rP;
		if(rC==0.0)
		{
			rP = ref.missing(p.first,options.markovOrder);
		} else {
			rP = refmarkov[p.first].prob;
		}

		

		double qCt = qC - (qN*qP);
		double rCt = rC - (rN*rP);
		//cout << p.first << " Q " << qCt << " from " << qC << ":" << qN << " " << qP << " " << (qN*qP) << endl;
		//cout << p.first << " R " << rCt << " from " << rC << ":" << rN << " " << rP << " " << (rN*rP) << endl;
		double dist = sqrt(qCt*qCt + rCt*rCt);

		D2S = D2S + (qCt*rCt / dist);
		sum1 = sum1 + (qCt*qCt / dist);
		sum2 = sum2 + (rCt*rCt / dist);
	}

	//cout << "rN,qN" << rN << " " << qN << endl;
	//cout << "rtot,qtot" << rtot << " " << qtot << endl;
	score = 0.5 * (1 - D2S/( sqrt(sum1)*sqrt(sum2) ));
	return score;
}

/*
double d2star(Seq ref, Seq query, ModifyStringOptions options)
{
	
}
*/
