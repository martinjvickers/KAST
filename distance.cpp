#include "distance.h"

double d2s(ModifyStringOptions options, map<string, bool> ourkmers, map<string, unsigned int> refcounts, map<string, double> refmarkov, map<string, unsigned int> querycounts, map<string, double> querymarkov)
{
        double score = 0.0;
        double D2S = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;

        int qN = 0;
        int rN = 0;

        for(pair<string, bool> p: ourkmers)
        {
                qN = qN + querycounts[p.first];
                rN = rN + refcounts[p.first];
        }

        for(pair<string, bool> p: ourkmers)
        {
                double qC = querycounts[p.first];
                double qP = querymarkov[p.first];
                double rC = refcounts[p.first];
                double rP = refmarkov[p.first];

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

double d2star(ModifyStringOptions options, map<string, bool> ourkmers, map<string, unsigned int> refcounts, map<string, double> refmarkov, map<string, unsigned int> querycounts, map<string, double> querymarkov)
{
        double D_2Star = 0.0;
        double tempQ = 0.0;
        double tempR = 0.0;

        int qN = 0;
        int rN = 0;

        for(pair<string, bool> p: ourkmers)
        {
                qN = qN + querycounts[p.first];
                rN = rN + refcounts[p.first];
        }

        for(pair<string, bool> p: ourkmers)
        {
                double qC = querycounts[p.first];
                double qP = querymarkov[p.first];
                double rC = refcounts[p.first];
                double rP = refmarkov[p.first];

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
}

double hao(ModifyStringOptions options, map<string, bool> ourkmers, map<string, unsigned int> refcounts, map<string, double> refmarkov, map<string, unsigned int> querycounts, map<string, double> querymarkov)
{
        int qN = 0;
        int rN = 0;

        for(pair<string, bool> p: ourkmers)
        {
                qN = qN + querycounts[p.first];
                rN = rN + refcounts[p.first];
        }

        double tempX = 0.0;
        double tempY = 0.0;
        double tempXY = 0.0;
	double d_Hao = 0.0;

        for(pair<string, bool> p: ourkmers)
        {
                double qC = querycounts[p.first];
                double qP = querymarkov[p.first];
                double rC = refcounts[p.first];
                double rP = refmarkov[p.first];

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
}

double euler(ModifyStringOptions options, map<string, unsigned int> refcounts, map<string, unsigned int> querycounts)
{
        double score = 0.0;
        int rN = 0;
        int qN = 0;

        for(pair<string, unsigned int> p: refcounts)
                rN = rN + refcounts[p.first];

        for(pair<string, unsigned int> p: querycounts)
                qN = qN + querycounts[p.first];

        //create a unified map
        map<string, unsigned int> ourkmers;
        ourkmers = refcounts;
        ourkmers.insert(querycounts.begin(),querycounts.end());

        for(pair<string, unsigned int> p: ourkmers)
        {
                double rF = refcounts[p.first] / (double)rN;
                double qF = querycounts[p.first] / (double)qN;
                score = score + (pow((rF - qF), 2));
        }

        return pow(score, 0.5);
}

double d2(ModifyStringOptions options, map<string, unsigned int> refcounts, map<string, unsigned int> querycounts)
{
        double sumqCrC = 0.0;
        double sumqC2 = 0.0;
        double sumrC2 = 0.0;

        //create a unified map
        map<string, unsigned int> ourkmers;
        ourkmers = refcounts;
        ourkmers.insert(querycounts.begin(),querycounts.end());

        for(pair<string, unsigned int> p: ourkmers)
        {
                sumqCrC = sumqCrC + (refcounts[p.first] * querycounts[p.first]);
                sumqC2 = sumqC2 + (querycounts[p.first] * querycounts[p.first]);
                sumrC2 = sumrC2 + (refcounts[p.first] * refcounts[p.first]);
        }

        double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}

double manhattan(ModifyStringOptions options, map<string, unsigned int> refcounts, map<string, unsigned int> querycounts)
{
        double score = 0.0;
        int rN = 0;
        int qN = 0;

        for(pair<string, unsigned int> p: refcounts)
                rN = rN + refcounts[p.first];

        for(pair<string, unsigned int> p: querycounts)
                qN = qN + querycounts[p.first];

        //create a unified map
        map<string, unsigned int> ourkmers;
        ourkmers = refcounts;
        ourkmers.insert(querycounts.begin(), querycounts.end());

        for(pair<string, unsigned int> p: ourkmers)
        {
                double rF = refcounts[p.first] / (double)rN;
                double qF = querycounts[p.first] / (double)qN;
                score = score + (abs(rF - qF));
        }

	return score;
}

double chebyshev(ModifyStringOptions options, map<string, unsigned int> refcounts, map<string, unsigned int> querycounts)
{
        double score = 0.0;
        int rN = 0;
        int qN = 0;
	double temp = 0.0;

        for(pair<string, unsigned int> p: refcounts)
                rN = rN + refcounts[p.first];

        for(pair<string, unsigned int> p: querycounts)
                qN = qN + querycounts[p.first];

        //create a unified map
        map<string, unsigned int> ourkmers;
        ourkmers = refcounts;
        ourkmers.insert(querycounts.begin(),querycounts.end());

        for(pair<string, unsigned int> p: ourkmers)
        {
                double rF = refcounts[p.first] / (double)rN;
                double qF = querycounts[p.first] / (double)qN;
                temp = abs(rF - qF);
                if(temp > score)
                        score = temp;
        }

        return score;
}
