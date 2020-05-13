#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define e 2.718281828459045
#define pi 3.141592653589793
#define alpha 418.982887

/* Parametrii des folositi */
#define D 50 /* Numar de parametrii care trebuie optimizati*/
#define runTimes 100/* De cate ori e rulat algoritmul*/

/* Functii pentru testat*/
double Rosenbrock(double sol[D]);
double Griewank(double sol[D]);
double Rastrigin(double sol[D]);
double Ackley(double sol[D]);
double Schwefel(double sol[D]); /* nu merge */

typedef double (*FunctionCallback)(double sol[D]);

FunctionCallback chosenFunction = &Rosenbrock;

#define colonyDimension 40 /* Dimensiunea coloniei*/
#define foodNumber colonyDimension / 2 /* Numarul resurselor de mancare = jumatate din dimensiunea coloniei*/
#define limit 100 /* Limita pentru counter-ul de imbunatatire*/
#define maxCycle 3000 /* Nr de cicluri pentru cautarea hranei*/
#define lowerLimit -5.12 /* Limita inferioara */
#define upperLimit 5.12 /* Limita superioara*/

double foods[foodNumber][D]; /* Populatia mancarii*/
double foodsF[foodNumber]; /* Tine valoriile functiei mancarii*/
double fitness[foodNumber]; /* Tine valorile fitness-ului pentru fiecare mancare*/
double trial[foodNumber]; /* Retine valorile unde solutiile nu pot fi imbunatatite (counter-ul de imbunatatire)*/
double prob[foodNumber]; /* Contine probabilitatile ca o sursa de mancare sa fie aleasa*/
double solution[D]; /* Noua solutie produsa v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj}), j - parametru random, k - solutie random, diferita de i*/
double solF; /* Valoarea functiei a noii solutii*/
double fitnessSol; /* Valoarea fitness-ului noii solutii*/
int neighbour; /* Corespunde k-ului din ecuatia v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
int changedParam; /* Corespunde lui j*/
double minGlobal; /* Optimum*/
double paramGlobal[D]; /* Parametrii solutiei optium*/
double minsGlobal[runTimes]; /* Tine minGlobal pentru mai multe run-uri*/
double rnd; /* numar random intre [0, 1)*/

double CalculateFitness(double fun); /* Functia pentru fitness*/
void MemorizeBestSource(); /* Retine cea mai buna sursa de mancare*/
void init(int index); /* Initializeaza valorile initiale*/
void initial(); /* Initializeaza sursele de mancare*/
void SendEmployedBees(); /* Trimite albinele lucratoare */
void CalculateProbabilities(); /* Calculeaza probabilitatile */
void SendOnlookerBees();
void SendScoutBees();

int main()
{
    int iter, run, j;
    double mean;
    mean = 0;

    for (run = 0; run < runTimes; run++) {
        initial();
        MemorizeBestSource();
        for (iter = 0; iter < maxCycle; iter++) {
            SendEmployedBees();
            CalculateProbabilities();
            SendOnlookerBees();
            MemorizeBestSource();
            SendScoutBees();
        }
        for (j = 0; j < D; j++) {
            printf("GlobalParam[%d]: %f\n", j + 1, paramGlobal[j]);
        }
        printf("%d. run: %e \n", run + 1, minGlobal);
        minsGlobal[run] = minGlobal;
        mean = mean + minGlobal;
    }

    mean = mean / runTimes;
    printf("Means of %d runs: %e\n", runTimes, mean);
}

double Rosenbrock(double sol[D]) {
    int j;
    double top = 0;
    for (j = 0; j < D - 1; j++) {
        top = top + 100 * pow((sol[j + 1] - pow((sol[j]), (double)2)), (double)2) + pow((sol[j] - 1), (double)2);
    }
    return top;
}

double Griewank(double sol[D]) {
    int j;
    double top1, top2, top;
    top = 0;
    top1 = 0;
    top2 = 1;
    for (j = 0; j < D; j++) {
        top1 = top1 + pow((sol[j]), (double)2);
        top2 = top2 * cos((((sol[j]) / sqrt((double)(j + 1))) * M_PI) / 180);
    }
    top = (1 / (double)4000) * top1 - top2 + 1;
    return top;
}

double Rastrigin(double sol[D]) {
    int j;
    double top = 0;

    for (j = 0; j < D; j++) {
        top = top + (pow(sol[j], (double)2) - 10 * cos(2 * M_PI * sol[j]) + 10);
    }
    return top;
}

double Ackley(double sol[D]) {
    int j;
    double top = 0;
    for (j = 0 ; j < D ; j++) {
        top = top + -20.0*exp(-0.2*sqrt(pow(sol[j], 2) / D )) - exp(cos(2 * pi * sol[D]) / D ) + 20 + e;
    }
    return top;
}

/* Nu merge */
double Schwefel(double sol[D]) {
    int j;
    double top = 0;
    for (j = 0 ; j < D ; j++) {
        top +=  sol[j] * sin(sqrt(fabs(sol[j])));
    }
    return (-1) * top + alpha * D;
}

double CalculateFitness(double fun)
{
    double result = 0;
    if (fun >= 0) {
        result = 1 / (fun + 1);
    }
    else {
        result = 1 + fabs(fun);
    }
    return result;
}

void MemorizeBestSource()
{
    int i, j;

    for (i = 0; i < foodNumber; i++) {
        if (foodsF[i] < minGlobal) {
            minGlobal = foodsF[i];
            for (j = 0; j < D; j++)
                paramGlobal[j] = foods[i][j];
        }
    }
}

void init(int index)
{
    int j;
    for (j = 0; j < D; j++) {
        rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        foods[index][j] = rnd * (upperLimit - lowerLimit) + lowerLimit;
        solution[j] = foods[index][j];
    }
    foodsF[index] = chosenFunction(solution);
    fitness[index] = CalculateFitness(foodsF[index]);
    trial[index] = 0;
}

void initial()
{
    int i;
    for (i = 0; i < foodNumber; i++) {
        init(i);
    }
    minGlobal = foodsF[0];
    for (i = 0; i < D; i++)
        paramGlobal[i] = foods[0][i];
}

void SendEmployedBees()
{
    int i, j;
    /* Etapa albinelor lucratoare*/
    for (i = 0; i < foodNumber; i++) {
        /* Parametrul care trebuie schimbat este ales aleator*/
        rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        changedParam = (int)(rnd * D);

        /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
        rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        neighbour = (int)(rnd * foodNumber);

        /* Solutia aleasa random trebuie sa fie diferita de i*/
        while (neighbour == i) {
            rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            neighbour = (int)(rnd * foodNumber);
        }
        for (j = 0; j < D; j++) {
            solution[j] = foods[i][j];
        }

        /* v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj})*/
        rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        solution[changedParam] = foods[i][changedParam] + (foods[i][changedParam] - foods[neighbour][changedParam]) * (rnd - 0.5) * 2;

        /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
        if (solution[changedParam] < lowerLimit)
            solution[changedParam] = lowerLimit;
        if (solution[changedParam] > upperLimit)
            solution[changedParam] = upperLimit;
        solF = chosenFunction(solution);
        fitnessSol = CalculateFitness(solF);

        if (fitnessSol > fitness[i]) {
            /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
            trial[i] = 0;
            for (j = 0; j < D; j++)
                foods[i][j] = solution[j];
            foodsF[i] = solF;
            fitness[i] = fitnessSol;
        }
        else {
            /* Daca solutia nu poate fi imbunatatita, creste counter-ul de imbunatatire*/
            trial[i] = trial[i] + 1;
        }
    }
}

/* O sursa de mancare este aleasa cu o probabilitate propotionala cu calitatea ei*/
/* Probabilitatile sunt calculate folosind valoarea de fitness si normalizata prin impartirea acesteia la valoarea maxima a fitness-ului*/
void CalculateProbabilities()
{
    int i;
    double maxfit;
    maxfit = fitness[0];
    for (i = 1; i < foodNumber; i++) {
        if (fitness[i] > maxfit)
            maxfit = fitness[i];
    }

    for (i = 0; i < foodNumber; i++) {
        prob[i] = (0.9 * (fitness[i] / maxfit)) + 0.1;
    }
}

void SendOnlookerBees()
{
    int i, j, t;
    i = 0;
    t = 0;
    /* Etapa albinelor care aleg sursa de mancare*/
    while (t < foodNumber) {

        rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        if (rnd < prob[i]) /* Alege o sursa de mancare in functie de probabilitatea ei de a fi aleasa*/
        {
            t++;

            /* Parametrul care trebuie schimbat este ales aleator*/
            rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            changedParam = (int)(rnd * D);

            /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
            rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            neighbour = (int)(rnd * foodNumber);

            /* Solutia aleasa random trebuie sa fie diferita de i*/
            while (neighbour == i) {
                rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
                neighbour = (int)(rnd * foodNumber);
            }
            for (j = 0; j < D; j++)
                solution[j] = foods[i][j];

            /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
            rnd = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            solution[changedParam] = foods[i][changedParam] + (foods[i][changedParam] - foods[neighbour][changedParam]) * (rnd - 0.5) * 2;

            /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
            if (solution[changedParam] < lowerLimit)
                solution[changedParam] = lowerLimit;
            if (solution[changedParam] > upperLimit)
                solution[changedParam] = upperLimit;
            solF = chosenFunction(solution);
            fitnessSol = CalculateFitness(solF);

            if (fitnessSol > fitness[i]) {
                /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
                trial[i] = 0;
                for (j = 0; j < D; j++)
                    foods[i][j] = solution[j];
                foodsF[i] = solF;
                fitness[i] = fitnessSol;
            }
            else {
                /* Daca solutia nu poate fi imbunatatita, creste counter-ul de imbunatatire*/
                trial[i] = trial[i] + 1;
            }
        }
        i++;
        if (i == foodNumber)
            i = 0;
    }
}

/* Detecteaza sursele pentru care counter-ul de imbunatatire depaseste limita si o abandoneaza*/
void SendScoutBees()
{
    int maxtrialindex, i;
    maxtrialindex = 0;
    for (i = 1 ; i < foodNumber; i++) {
        if (trial[i] > trial[maxtrialindex])
            maxtrialindex = i;
    }
    if (trial[maxtrialindex] >= limit) {
        init(maxtrialindex);
    }
}
