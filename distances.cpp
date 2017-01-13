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

        //create union of the two maps
        unordered_map<string, markov_dat> ourkmers;
        ourkmers = refmarkov;
        ourkmers.insert(qrymarkov.begin(),qrymarkov.end());

        for(pair<string, markov_dat> p: ourkmers)
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

double d2sopt(Seq ref, Seq query, ModifyStringOptions options, map<string, bool> ourkmers)
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

