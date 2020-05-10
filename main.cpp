#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

/* Constante*/
#define e 2.718281828459045
#define pi 3.141592653589793
#define alpha 418.982887

/* Parametrii de control*/
#define NP 40 /* Dimensiunea coloniei*/
#define FoodNumber NP / 2 /* Numarul resurselor de mancare = jumatate din dimensiunea coloniei*/
#define limit 100 /* Limita pentru counter-ul de imbunatatire*/
#define maxCycle 3000 /* Nr de cicluri pentru cautarea hranei*/

/* Constante specifice problemei*/
#define D 100 /* Numar de parametrii care trebuie optimizati*/
#define lb -5.12 /* Limita inferioara */
#define ub 5.12 /* Limita superioara*/

#define runtime 30 /* De cate ori e rulat algoritmul*/

double Foods[FoodNumber][D]; /* Populatia mancarii*/
double f[FoodNumber]; /* Tine valoriile functiei mancarii*/
double fitness[FoodNumber]; /* Tine valorile fitness-ului pentru fiecare mancare*/
double trial[FoodNumber]; /* Retine valorile unde solutiile nu pot fi imbunatatite (counter-ul de imbunatatire)*/
double prob[FoodNumber]; /* Contine probabilitatile ca o sursa de mancare sa fie aleasa*/
double solution[D]; /* Noua solutie produsa v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj}), j - parametru random, k - solutie random, diferita de i*/
double ObjValSol; /* Valoarea functiei a noii solutii*/
double FitnessSol; /* Valoarea fitness-ului noii solutii*/
int neighbour; /* Corespunde k-ului din ecuatia v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
int param2change; /* Corespunde lui j*/
double GlobalMin; /* Optimum*/
double GlobalParams[D]; /* Parametrii solutiei optium*/
double GlobalMins[runtime]; /* Tine GlobalMin pentru mai multe run-uri*/
double r; /* numar random intre [0, 1)*/

/* Settup pentru functia folosita*/
typedef double (*FunctionCallback)(double sol[D]);

/* Functii pentru testat*/
double Rosenbrock(double sol[D]);
double Griewank(double sol[D]);
double Rastrigin(double sol[D]);
double Ackley(double sol[D]);
double Schwefel(double sol[D]);

/* Functia folosita*/
FunctionCallback function = &Ackley;

/* Functia pentru fitness*/
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

/* Retine cea mai buna sursa de mancare*/
void MemorizeBestSource()
{
    int i, j;

    for (i = 0; i < FoodNumber; i++) {
        if (f[i] < GlobalMin) {
            GlobalMin = f[i];
            for (j = 0; j < D; j++)
                GlobalParams[j] = Foods[i][j];
        }
    }
}

/* Initializeaza valorile initiale*/
void init(int index)
{
    int j;
    for (j = 0; j < D; j++) {
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        Foods[index][j] = r * (ub - lb) + lb;
        solution[j] = Foods[index][j];
    }
    f[index] = function(solution);
    fitness[index] = CalculateFitness(f[index]);
    trial[index] = 0;
}

/* Initializeaza sursele de mancare*/
void initial()
{
    int i;
    for (i = 0; i < FoodNumber; i++) {
        init(i);
    }
    GlobalMin = f[0];
    for (i = 0; i < D; i++)
        GlobalParams[i] = Foods[0][i];
}

void SendEmployedBees()
{
    int i, j;
    /* Etapa albinelor lucratoare*/
    for (i = 0; i < FoodNumber; i++) {
        /* Parametrul care trebuie schimbat este ales aleator*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        param2change = (int)(r * D);

        /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        neighbour = (int)(r * FoodNumber);

        /* Solutia aleasa random trebuie sa fie diferita de i*/
        while (neighbour == i) {
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            neighbour = (int)(r * FoodNumber);
        }
        for (j = 0; j < D; j++) {
            solution[j] = Foods[i][j];
        }

        /* v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj})*/
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        solution[param2change] = Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbour][param2change]) * (r - 0.5) * 2;

        /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
        if (solution[param2change] < lb)
            solution[param2change] = lb;
        if (solution[param2change] > ub)
            solution[param2change] = ub;
        ObjValSol = function(solution);
        FitnessSol = CalculateFitness(ObjValSol);

        if (FitnessSol > fitness[i]) {
            /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
            trial[i] = 0;
            for (j = 0; j < D; j++)
                Foods[i][j] = solution[j];
            f[i] = ObjValSol;
            fitness[i] = FitnessSol;
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
    for (i = 1; i < FoodNumber; i++) {
        if (fitness[i] > maxfit)
            maxfit = fitness[i];
    }

    for (i = 0; i < FoodNumber; i++) {
        prob[i] = (0.9 * (fitness[i] / maxfit)) + 0.1;
    }
}

void SendOnlookerBees()
{
    int i, j, t;
    i = 0;
    t = 0;
    /* Etapa albinelor care aleg sursa de mancare*/
    while (t < FoodNumber) {

        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        if (r < prob[i]) /* Alege o sursa de mancare in functie de probabilitatea ei de a fi aleasa*/
        {
            t++;

            /* Parametrul care trebuie schimbat este ales aleator*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            param2change = (int)(r * D);

            /* O solutie aleasa random creaza o solutie derivata a solutiei i*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            neighbour = (int)(r * FoodNumber);

            /* Solutia aleasa random trebuie sa fie diferita de i*/
            while (neighbour == i) {
                r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
                neighbour = (int)(r * FoodNumber);
            }
            for (j = 0; j < D; j++)
                solution[j] = Foods[i][j];

            /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
            r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
            solution[param2change] = Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbour][param2change]) * (r - 0.5) * 2;

            /* Daca valoarea este inafara limitelor se seteaza ca valoare limita*/
            if (solution[param2change] < lb)
                solution[param2change] = lb;
            if (solution[param2change] > ub)
                solution[param2change] = ub;
            ObjValSol = function(solution);
            FitnessSol = CalculateFitness(ObjValSol);

            if (FitnessSol > fitness[i]) {
                /* Daca solutia derivata este mai buna decat cea curenta, noua solutie devine cea derivata si reseteaza counter-ul de imbunatatire*/
                trial[i] = 0;
                for (j = 0; j < D; j++)
                    Foods[i][j] = solution[j];
                f[i] = ObjValSol;
                fitness[i] = FitnessSol;
            }
            else {
                /* Daca solutia nu poate fi imbunatatita, creste counter-ul de imbunatatire*/
                trial[i] = trial[i] + 1;
            }
        }
        i++;
        if (i == FoodNumber)
            i = 0;
    }
}

/* Detecteaza sursele pentru care counter-ul de imbunatatire depaseste limita si o abandoneaza*/
void SendScoutBees()
{
    int maxtrialindex, i;
    maxtrialindex = 0;
    for (i = 1 ; i < FoodNumber; i++) {
        if (trial[i] > trial[maxtrialindex])
            maxtrialindex = i;
    }
    if (trial[maxtrialindex] >= limit) {
        init(maxtrialindex);
    }
}

int main()
{
    int iter, run, j;
    double mean;
    mean = 0;

    for (run = 0; run < runtime; run++) {
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
            printf("GlobalParam[%d]: %f\n", j + 1, GlobalParams[j]);
        }
        printf("%d. run: %e \n", run + 1, GlobalMin);
        GlobalMins[run] = GlobalMin;
        mean = mean + GlobalMin;
    }

    mean = mean / runtime;
    printf("Means of %d runs: %e\n", runtime, mean);
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

/* Nu merge*/
double Schwefel(double sol[D]) {
    int j;
    double top = 0;
    for (j = 0 ; j < D ; j++) {
        top = top - sol[j] * sin(sqrt(fabs(sol[j])));
    }
    return top + alpha * D;
}