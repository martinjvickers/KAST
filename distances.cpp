#include "common.h"
#include "distances.h"

double d2(long long int ref[], long long int qry[], int nokmers)
{
        double sumqCrC = 0.0;
        double sumqC2 = 0.0;
        double sumrC2 = 0.0;

        for(int i = 0; i < nokmers; i++){
                sumqCrC = sumqCrC + (qry[i] * ref[i]);
                sumqC2 = sumqC2 + (qry[i] * qry[i]);
                sumrC2 = sumrC2 + (ref[i] * ref[i]);
        }

        double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
        return 0.5*(1-score);
}

double euler(long long int ref[], long long int qry[], int nokmers)
{
        double score = 0.0;
        double rN = 0.0;
        double qN = 0.0;

        for(int i = 0; i < nokmers; i++){
                rN = rN + ref[i];
                qN = qN + qry[i];
        }

        for(int i = 0; i < nokmers; i++){
                double rF = ref[i] / rN;
                double qF = qry[i] / qN;
                score = score + (pow((rF - qF), 2));
        }

        return pow(score, 0.5);
}

double d2s(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[])
{
        double score = 0.0;
        double D2S = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double rN = 0.0;
        double qN = 0.0;
        double rtot = 0.0;
        double qtot = 0.0;

        for(int i = 0; i < nokmers; i++)
        {
                qN = qN + qryCounts[i];
                qtot = qtot + qryProbs[i];
                rN = rN + refCounts[i];
                rtot = rtot + refProbs[i];
        }

        for(int i = 0; i < nokmers; i++)
        {
                double qC = qryCounts[i];
                double rC = refCounts[i];
                double qP = qryProbs[i];
                double rP = refProbs[i];

                double qCt = qC - (qN * qP);
                double rCt = rC - (rN * rP);
                double dist = sqrt(( qCt * qCt )+( rCt * rCt ));

                if(dist == 0.0)
                {
                        std::cout << "Div by zero" << std::endl;
                }

                D2S = D2S + ( ( qCt * rCt ) / dist );
                sum1 = sum1 + ( ( qCt * qCt ) / dist );
                sum2 = sum2 + ( ( rCt * rCt ) / dist );
        }

        score = 0.5 * (1 - ( (D2S) / (sqrt(sum1)*sqrt(sum2)) ) );

        return score;

}

double d2star(long long int refCounts[], long long int qryCounts[], int nokmers, double refProbs[], double qryProbs[])
{
        double score = 0.0;
        double D2Star = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double qN = 0.0;
        double rN = 0.0;

        for(int i = 0; i < nokmers; i++)
        {
                qN = qN + qryCounts[i];
                rN = rN + refCounts[i];
        }

        for(int i = 0; i < nokmers; i++)
        {
                double qC = qryCounts[i];
                double rC = refCounts[i];
                double qP = qryProbs[i];
                double rP = refProbs[i];

                double q_np = qN * qP;
                double r_np = rN * rP;
                double qCt = qC - q_np;
                double rCt = rC - r_np;

                double dist = sqrt(q_np) * sqrt(r_np);

                if(dist == 0.0)
                {
                        std::cout << "Div by zero" << std::endl;
                }

                D2Star = D2Star + ( ( qCt * rCt ) / dist );
                sum1 = sum1 + ( ( qCt * qCt ) / q_np );
                sum2 = sum2 + ( ( rCt * rCt ) / r_np );
        }
        score = 0.5 * (1 - ( (D2Star) / (sqrt(sum1)*sqrt(sum2)) ) );
        return score;
}
