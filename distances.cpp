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

double manhattan(Seq ref, Seq query, ModifyStringOptions options)
{
	double score = 0.0;

	//extract our maps
	unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
	unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

	//total number of counts for query and ref
	int rN = ref.getTotalCounts(options.klen);
	int qN = query.getTotalCounts(options.klen);

	//create a unified map
        unordered_map<string, long long int> ourkmers;
        ourkmers = refkmers;
        ourkmers.insert(querykmers.begin(), querykmers.end());

	for(pair<string, long long int> p: ourkmers)
        {
		double rF = refkmers[p.first] / rN;
		double qF = querykmers[p.first] / qN;
		score = score + (abs(rF - qF));
	}

	return pow(score, 0.5);

}

double chebyshev(Seq ref, Seq query, ModifyStringOptions options)
{
	double score = 0.0;

	//extract our maps
        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

	//total number of counts for query and ref
        int rN = ref.getTotalCounts(options.klen);
        int qN = query.getTotalCounts(options.klen);

	//create a unified map
        unordered_map<string, long long int> ourkmers;
        ourkmers = refkmers;
        ourkmers.insert(querykmers.begin(), querykmers.end());

	for(pair<string, long long int> p: ourkmers)
        {
		double rF = refkmers[p.first] / rN;
                double qF = querykmers[p.first] / qN;
		double temp = abs(rF - qF);
		if(temp > score)
			score = temp;
	}
	return score;
}

double d2s(Seq ref, Seq query, ModifyStringOptions options, map<string, bool> ourkmers)
{
        int rN = ref.getTotalCounts(options.klen);
        int qN = query.getTotalCounts(options.klen);
        double score = 0.0;
        double D2S = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;

        //get my markov maps
        unordered_map<string,markov_dat> qrymarkov = query.getMarkov(options.klen, options.markovOrder);
        unordered_map<string,markov_dat> refmarkov = ref.getMarkov(options.klen, options.markovOrder);

        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);

        for(pair<string, bool> p: ourkmers)
        {
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
                double dist = sqrt(qCt*qCt + rCt*rCt);

                D2S = D2S + (qCt*rCt / dist);
                sum1 = sum1 + (qCt*qCt / dist);
                sum2 = sum2 + (rCt*rCt / dist);
        }
        score = 0.5 * (1 - D2S/( sqrt(sum1)*sqrt(sum2) ));
        return score;
}

double d2star(Seq ref, Seq query, ModifyStringOptions options, map<string, bool> ourkmers)
{
	double nX = query.getTotalCounts(options.klen);
	double nY = ref.getTotalCounts(options.klen);

        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);
        unordered_map<string,markov_dat> qrymarkov = query.getMarkov(options.klen, options.markovOrder);
        unordered_map<string,markov_dat> refmarkov = ref.getMarkov(options.klen, options.markovOrder);

	double D_2Star = 0.0;
        double tempX = 0.0;
        double tempY = 0.0;

	for(pair<string, bool> p: ourkmers)
        {
		double cXi = querykmers[p.first]; //get the counts
		double pXi;
		if(cXi == 0.0)
		{
			pXi = query.missing(p.first, options.markovOrder);
		} else {
			pXi = qrymarkov[p.first].prob;
		}

		double cYi = refkmers[p.first];
		double pYi;
		if(refkmers.find(p.first) != refkmers.end())
		{
			pYi = ref.missing(p.first, options.markovOrder);
		} else {
			pYi = refmarkov[p.first].prob;
		}

		double temp1 = nX * pXi;
		double temp2 = nY * pYi;
		double cXi_bar = cXi - temp1;
		double cYi_bar = cYi - temp2;
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
	
		D_2Star += cXi_bar * cYi_bar / temp3;
		tempX += cXi_bar * cXi_bar / temp1;
		tempY += cYi_bar * cYi_bar / temp2;
	}

	double temp = D_2Star / (sqrt(tempX) * sqrt(tempY));
	return 0.5 * (1 - temp);
}

double hao(Seq ref, Seq query, ModifyStringOptions options, map<string, bool> ourkmers)
{
	double nX = query.getTotalCounts(options.klen);
	double nY = ref.getTotalCounts(options.klen);

        unordered_map<string, long long int> querykmers = query.getCounts(options.klen);
        unordered_map<string, long long int> refkmers = ref.getCounts(options.klen);
        unordered_map<string,markov_dat> qrymarkov = query.getMarkov(options.klen, options.markovOrder);
        unordered_map<string,markov_dat> refmarkov = ref.getMarkov(options.klen, options.markovOrder);

	double tempX = 0.0;
	double tempY = 0.0;
	double tempXY = 0.0;

	for(pair<string, bool> p: ourkmers) //for each kmer
        {

		double cXi = querykmers[p.first]; //get the counts
		double pXi;
                if(cXi == 0.0)
                {
                        pXi = query.missing(p.first, options.markovOrder);
                } else {
                        pXi = qrymarkov[p.first].prob;
                }

                double cYi = refkmers[p.first];
                double pYi;
                if(refkmers.find(p.first) != refkmers.end())
                {
                        pYi = ref.missing(p.first, options.markovOrder);
                } else {
                        pYi = refmarkov[p.first].prob;
                }
	
		double fXi = cXi / nX;
		double fYi = cYi / nY;
		double temp1, temp2;

		if(pXi == 0)
			temp1 = -1;
		else
			temp1 = (fXi / pXi) - 1.0;

		if(pYi == 0)
			temp2 = -1;
		else
			temp2 = (fYi / pYi) - 1.0;

		tempXY += temp1 * temp2;
		tempX += temp1 * temp1;
		tempY += temp2 * temp2;
	}

	double temp = tempXY / (sqrt(tempX) * sqrt(tempY));
	return (1 - temp) / 2;

}
