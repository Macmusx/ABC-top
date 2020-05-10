#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Constante*/
#define e 2.718281828459045
#define pi 3.141592653589793
#define alpha 418.982887

/* Parametrii de control*/
#define columnDimension 40 /* Dimensiunea coloniei*/
#define numberOfFoodResources columnDimension / 2 /* Numarul resurselor de mancare = jumatate din dimensiunea coloniei*/
#define limit 100 /* Limita pentru counter-ul de imbunatatire*/
#define maxNumberOfCycles 3000 /* Numarul de cicluri pentru cautarea hranei*/

/* Constante specifice problemei*/
#define noOfOptmimizedParams 50 /* Numar de parametrii care trebuie optimizati*/
#define lowLimit -5.12 /* Limita inferioara */
#define upperLimit 5.12 /* Limita superioara*/

#define numberOfRuns 10 /* De cate ori este rulat algoritmul*/

double Foods[numberOfFoodResources][noOfOptmimizedParams]; /* Populatia mancarii*/
double f[numberOfFoodResources]; /* Stocheaza valorile functiei mancarii*/
double fitness[numberOfFoodResources]; /* Tine valorile fitness-ului pentru fiecare mancare*/
double improvementFactor[numberOfFoodResources]; /* Retine valorile unde solutiile nu pot fi imbunatatite (counter-ul de imbunatatire)*/
double probForFoodResToBeChoosen[numberOfFoodResources]; /* Contine probabilitatile ca o sursa de mancare sa fie aleasa*/
double solution[noOfOptmimizedParams]; /* Noua solutie produsa v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj}), j - parametru random, k - solutie random, diferita de i*/
double functionNewVal; /* Valoarea functiei noii solutii*/
double fitnessVal; /* Valoarea fitness-ului noii solutii*/
int proximity; /* Corespunde k din ecuatia v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
int paramForChange; /* Parametrul care trebuie schimbat, corespunde lui j*/
double GlobalMin; /* Minimul global - valoarea optima*/
double GlobalParams[noOfOptmimizedParams]; /* Parametrii solutiei minimului global ( a valorii optime )*/
double GlobalMins[numberOfRuns]; /* Stocheaza valorile pentru GlobalMin pentru fiecare rulare*/
double r; /* numar random intre [0, 1)*/

/* Settup pentru functia folosita*/
typedef double (*FunctionCallback)(double sol[noOfOptmimizedParams]);

/* Functii pentru testat*/
double Rosenbrock(double sol[noOfOptmimizedParams]);
double Griewank(double sol[noOfOptmimizedParams]);
double Rastrigin(double sol[noOfOptmimizedParams]);
double Ackley(double sol[noOfOptmimizedParams]);

/* Functia folosita*/
FunctionCallback function = &Ackley;

/* Functia pentru fitness*/
double CalculateFitness(double index)
{
    double result = 0;
    if (index >= 0) {
        result = 1 / (index + 1);
    }
    else {
        result = 1 + fabs(index);
    }
    return result;
}

/* Retine cea mai buna sursa de mancare*/
void MemorizeBestSource()
{
    int i, j;

    for (i = 0; i < numberOfFoodResources; i++) {
        if (f[i] < GlobalMin) {
            GlobalMin = f[i];
            for (j = 0; j < noOfOptmimizedParams; j++)
                GlobalParams[j] = Foods[i][j];
        }
    }
}

/* Initializeaza valorile initiale*/
void init(int index)
{
    int j;
    for (j = 0; j < noOfOptmimizedParams; j++) {
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        Foods[index][j] = r * (upperLimit - lowLimit) + lowLimit;
        solution[j] = Foods[index][j];
    }
    f[index] = function(solution);
    fitness[index] = CalculateFitness(f[index]);
    improvementFactor[index] = 0;
}

/* Initializeaza sursele de mancare*/
void initial()
{
    int i;
    for (i = 0; i < numberOfFoodResources; i++) {
        init(i);
    }
    GlobalMin = f[0];
    for (i = 0; i < noOfOptmimizedParams; i++)
        GlobalParams[i] = Foods[0][i];
}

void SendEmployedBees()
{
    int i, j;
    /* Etapa albinelor lucratoare*/
    for (i = 0; i < numberOfFoodResources; i++) {
        /* Parametrul care trebuie schimbat este ales aleator*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        paramForChange = (int)(r * noOfOptmimizedParams);

        /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        proximity = (int)(r * numberOfFoodResources);

        /* Solutia aleasa random trebuie sa fie diferita de i*/
        while (proximity == i) {
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            proximity = (int)(r * numberOfFoodResources);
        }
        for (j = 0; j < noOfOptmimizedParams; j++) {
            solution[j] = Foods[i][j];
        }

        /* v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj})*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        solution[paramForChange] = Foods[i][paramForChange] + (Foods[i][paramForChange] - Foods[proximity][paramForChange]) * (r - 0.5) * 2;

        /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
        if (solution[paramForChange] < lowLimit)
            solution[paramForChange] = lowLimit;
        if (solution[paramForChange] > upperLimit)
            solution[paramForChange] = upperLimit;
        functionNewVal = function(solution);
        fitnessVal = CalculateFitness(functionNewVal);

        if (fitnessVal > fitness[i]) {
            /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
            improvementFactor[i] = 0;
            for (j = 0; j < noOfOptmimizedParams; j++)
                Foods[i][j] = solution[j];
            f[i] = functionNewVal;
            fitness[i] = fitnessVal;
        }
        else {
            /* Daca solutia nu poate fi imbunatatita, creste counter-ul de imbunatatire*/
            improvementFactor[i] = improvementFactor[i] + 1;
        }
    }
}

/* O sursa de mancare este aleasa cu o probabilitate propotionala cu calitatea ei*/
/* Probabilitatile sunt calculate folosind valoarea de fitness si normalizata prin impartirea acesteia la valoarea maxima a fitness-ului*/
void CalculateProbabilities()
{
    int i;
    double maxiumumFitness;
    maxiumumFitness = fitness[0];
    for (i = 1; i < numberOfFoodResources; i++) {
        if (fitness[i] > maxiumumFitness)
            maxiumumFitness = fitness[i];
    }

    for (i = 0; i < numberOfFoodResources; i++) {
        probForFoodResToBeChoosen[i] = (0.9 * (fitness[i] / maxiumumFitness)) + 0.1;
    }
}

void SendOnlookerBees()
{
    int i, j, t;
    i = 0;
    t = 0;
    /* Etapa albinelor cautatoare*/
    while (t < numberOfFoodResources) {

        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        if (r < probForFoodResToBeChoosen[i]) /* Alege o sursa de mancare in functie de probabilitatea acesteia de a fi aleasa*/
        {
            t++;

            /* Parametrul care trebuie schimbat este ales aleator*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            paramForChange = (int)(r * noOfOptmimizedParams);

            /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            proximity = (int)(r * numberOfFoodResources);

            /* Solutia aleasa random trebuie sa fie diferita de i*/
            while (proximity == i) {
                r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
                proximity = (int)(r * numberOfFoodResources);
            }
            for (j = 0; j < noOfOptmimizedParams; j++)
                solution[j] = Foods[i][j];

            /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            solution[paramForChange] = Foods[i][paramForChange] + (Foods[i][paramForChange] - Foods[proximity][paramForChange]) * (r - 0.5) * 2;

            /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
            if (solution[paramForChange] < lowLimit)
                solution[paramForChange] = lowLimit;
            if (solution[paramForChange] > upperLimit)
                solution[paramForChange] = upperLimit;
            functionNewVal = function(solution);
            fitnessVal = CalculateFitness(functionNewVal);

            /*a greedy selection is applied between the current solution i and its mutant*/
            if (fitnessVal > fitness[i]) {
                /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
                improvementFactor[i] = 0;
                for (int j = 0; j < noOfOptmimizedParams; j++)
                    Foods[i][j] = solution[j];
                f[i] = functionNewVal;
                fitness[i] = fitnessVal;
            }
            else {
                /* Daca solutia nu poate fi imbunatatita, creste counter-ul de imbunatatire*/
                improvementFactor[i] = improvementFactor[i] + 1;
            }
        }
        i++;
        if (i == numberOfFoodResources)
            i = 0;
    }
}

/* Detecteaza sursele pentru care counter-ul de imbunatatire depaseste limita si o abandoneaza*/
void SendScoutBees()
{
    int maxImprovementFactorIndex, i;
    maxImprovementFactorIndex = 0;
    for (i = 1 ; i < numberOfFoodResources; i++) {
        if (improvementFactor[i] > improvementFactor[maxImprovementFactorIndex])
            maxImprovementFactorIndex = i;
    }
    if (improvementFactor[maxImprovementFactorIndex] >= limit) {
        init(maxImprovementFactorIndex);
    }
}

int main()
{
    double globalMinMean;
    globalMinMean = 0;
    srand(time(NULL));

    for (int run = 0; run < numberOfRuns; run++) {
        initial();
        MemorizeBestSource();
        for (int iterator = 0; iterator < maxNumberOfCycles; iterator++) {
            SendEmployedBees();
            CalculateProbabilities();
            SendOnlookerBees();
            MemorizeBestSource();
            SendScoutBees();
        }
        for (int j = 0; j < noOfOptmimizedParams; j++) {
           printf("GlobalParam[%d]: %f\n", j + 1, GlobalParams[j]);
        }
        printf("%d. run: %e \n", run + 1, GlobalMin);
        GlobalMins[run] = GlobalMin;
        globalMinMean = globalMinMean + GlobalMin;
    }

    globalMinMean = globalMinMean / numberOfRuns;
    printf("Means of %d runs: %e\n", numberOfRuns, globalMinMean);
}

double Rosenbrock(double sol[noOfOptmimizedParams]) {
    double top = 0;
    for (int j = 0; j < noOfOptmimizedParams - 1; j++) {
        top = top + 100 * pow((sol[j + 1] - pow((sol[j]), (double)2)), (double)2) + pow((sol[j] - 1), (double)2);
    }
    return top;
}

double Griewank(double sol[noOfOptmimizedParams]) {
    double top, top1, top2;
    top = 0;
    top1 = 0;
    top2 = 1;
    for (int j = 0; j < noOfOptmimizedParams; j++) {
        top1 = top1 + pow((sol[j]), (double)2);
        top2 = top2 * cos((((sol[j]) / sqrt((double)(j + 1))) * M_PI) / 180);
    }
    top = (1 / (double)4000) * top1 - top2 + 1;
    return top;
}

double Rastrigin(double sol[noOfOptmimizedParams]) {
    double top = 0;

    for (int j = 0; j < noOfOptmimizedParams; j++) {
        top = top + (pow(sol[j], (double)2) - 10 * cos(2 * M_PI * sol[j]) + 10);
    }
    return top;
}

double Ackley(double sol[noOfOptmimizedParams]) {
    double top = 0;
    for (int j = 0 ; j < noOfOptmimizedParams ; j++) {
        top = top + -20.0*exp(-0.2*sqrt(pow(sol[j], 2) / noOfOptmimizedParams )) - exp(cos(2 * pi * sol[noOfOptmimizedParams]) / noOfOptmimizedParams ) + 20 + e;
    }
    return top;
}