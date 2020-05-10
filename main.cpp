#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>

#define rand_01 ((float)rand() / (float)RAND_MAX)

const int numofdims = 30;
const int numofparticles = 50;

using namespace std;

void fitnessfunc(float X[numofparticles][numofdims], float fitnesses[numofparticles])
{
    memset(fitnesses, 0, sizeof (float) * numofparticles);
    for(int i = 0; i < numofparticles; i++)
    {
        for(int j = 0; j < numofdims; j++)
        {
            fitnesses[i] += (pow(X[i][j], 2));
        }
    }
}

float mean(float inputval[], int vallength)
{
    int addvalue = 0;
    for(int i = 0; i < vallength; i++)
    {
        addvalue += inputval[i];
    }
    return (float)(addvalue / vallength);
}

void PSO(int numofiterations, float c1, float c2,
         float Xmin[numofdims], float Xmax[numofdims], float initialpop[numofparticles][numofdims],
         float worsts[], float meanfits[], float bests[], float *gbestfit, float gbest[numofdims])
{
    float V[numofparticles][numofdims] = {0};
    float X[numofparticles][numofdims];
    float Vmax[numofdims];
    float Vmin[numofdims];
    float pbests[numofparticles][numofdims];
    float pbestfits[numofparticles];
    float fitnesses[numofparticles];
    float w;
    float minfit;
    int   minfitidx;

    memcpy(X, initialpop, sizeof(float) * numofparticles * numofdims);
    fitnessfunc(X, fitnesses);
    minfit = *min_element(fitnesses, fitnesses + numofparticles);
    minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
    *gbestfit = minfit;
    memcpy(gbest, X[minfitidx], sizeof(float) * numofdims);

    for(int i = 0; i < numofdims; i++)
    {
        Vmax[i] = 0.2 * (Xmax[i] - Xmin[i]);
        Vmin[i] = -Vmax[i];
    }

    for(int t = 0; t < 1000; t++)
    {
        w = 0.9 - 0.7 * t / numofiterations;

        for(int i = 0; i < numofparticles; i++)
        {
            if(fitnesses[i] < pbestfits[i])
            {
                pbestfits[i] = fitnesses[i];
                memcpy(pbests[i], X[i], sizeof(float) * numofdims);
            }
        }
        for(int i = 0; i < numofparticles; i++)
        {
            for(int j = 0; j < numofdims; j++)
            {
                V[i][j] = min(max((w * V[i][j] + rand_01 * c1 * (pbests[i][j] - X[i][j])
                                   + rand_01 * c2 * (gbest[j] - X[i][j])), Vmin[j]), Vmax[j]);
                X[i][j] = min(max((X[i][j] + V[i][j]), Xmin[j]), Xmax[j]);
            }
        }

        fitnessfunc(X, fitnesses);
        minfit = *min_element(fitnesses, fitnesses + numofparticles);
        minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
        if(minfit < *gbestfit)
        {
            *gbestfit = minfit;
            memcpy(gbest, X[minfitidx], sizeof(float) * numofdims);
        }

        worsts[t] = *max_element(fitnesses, fitnesses + numofparticles);
        bests[t] = *gbestfit;
        meanfits[t] = mean(fitnesses, numofparticles);
    }


}

int main()
{
    time_t t;
    srand((unsigned) time(&t));

    float xmin[30], xmax[30];
    float initpop[50][30];
    float worsts[1000], bests[1000];
    float meanfits[1000];
    float gbestfit;
    float gbest[30];
    for(int i = 0; i < 30; i++)
    {
        xmax[i] = 100;
        xmin[i] = -100;
    }
    for(int i = 0; i < 50; i++)
        for(int j = 0; j < 30; j++)
        {
            initpop[i][j] = rand() % (100 + 100 + 1) - 100;
        }

    PSO(1000, 2, 2, xmin, xmax, initpop, worsts, meanfits, bests, &gbestfit, gbest);

    cout<<"fitness: "<<gbestfit<<endl;
    return 0;
}