#include <iostream>
#include <ctime>
#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include <cstring>

// C++11 specific header
#include <random>

#include "main.h"
#include "stack.h"
#include "shunting_yard.h"
#include "rpn.h"

using namespace std;

// cross chromosomes by averaging values
#define CROSSOVER_AVERAGE 1
// cross chromosomes by substituting values
#define CROSSOVER_SUBSTITUTION 2

// size of population
static int population_size;
// generation count (number of iterations)
static int generation_count_limit;

// base chance for certain chromosome to be crossed with another
static double crossover_base_chance_pct = 0.30;
// base chance for certain chromosome to randomly change one element
static double mutation_base_chance_pct = 0.10;

// genome size (variable count in supplied expression)
static int genome_size;
// bottom limit of parameter
static double parameter_limit_min = -20.0;
// top limit of parameter
static double parameter_limit_max = 20.0;

// crossover method used (default is averaging)
static int crossover_method = CROSSOVER_AVERAGE;

// expression stack used for evaluation
c_stack* exp_stack = NULL;
// map of variables to be used
std::list<char> variable_map;
// evaluation map with values for variables
double* eval_map;
// right side of equation (constant)
double rightside;

// default structure of chromosome
typedef struct
{
    // array of genetic elements
    double* genome;
    // precalculated fitness
    double fitness;
} chromozome;

// returns random number in interval <0.0;1.0>
inline double frand()
{
    static std::uniform_real<double> unif(0, 1);
    static std::default_random_engine re;
    return unif(re);

    // for pre C++11 compilers
    //return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

// calculates fitness of supplied chromosome using prebuilt evaluation stack
double fitness_function(chromozome* chr)
{
    int i = 0;
    // update values in evaluation stack
    for (std::list<char>::iterator itr = variable_map.begin(); itr != variable_map.end(); ++itr)
        eval_map[*itr] = chr->genome[i++];

    // evaluate left side
    double leftside = rpn_evaluate_stack(exp_stack, eval_map);

    // the lower the distance, the higher the fitness
    // this is achieved using its multiplicative inverse
    return 1.0/fabs(leftside - rightside);
}

// recalculates fitness of all chromosomes in population
void recalculate_fitness(chromozome** population, int popsize)
{
    for (int i = 0; i < popsize; i++)
        population[i]->fitness = fitness_function(population[i]);
}

// selects random chromozome using its fitness (the more fit, the higher the chance to be selected)
chromozome* select_random_chromozome(chromozome** population, int popsize, chromozome* except)
{
    // at first, update fitness to count with actual data
    recalculate_fitness(population, popsize);
    // then calculate sum of fitnesses
    double fitsum = 0.0;
    for (int i = 0; i < popsize; i++)
    {
        // we skip the "except" chromosome
        if (population[i] == except)
            continue;
        fitsum += population[i]->fitness;
    }

    // random from 0 to sum_of_fitnesses
    double rnd = frand()*fitsum;

    // now iterate through chromosomes and return the one the random number hit
    fitsum = 0.0;
    for (int i = 0; i < popsize; i++)
    {
        // again skip the "except" chromosome
        if (population[i] == except)
            continue;

        fitsum += population[i]->fitness;
        if (rnd < fitsum)
            return population[i];
    }

    return NULL;
}

// proceeds crossover between two chromosomes
void crossover_single(chromozome* parent, chromozome* child)
{
    // get random position in chromosome, move by one from both sides to avoid pure copying
    int pos = 1 + (rand() % (genome_size - 2));

    // there is 50/50 chance of copying left/right side
    if (rand() % 2 == 0)
    {
        // copy right side of crossover point
        for (int i = pos; i < genome_size; i++)
        {
            if (crossover_method == CROSSOVER_AVERAGE)
                child->genome[i] = (child->genome[i] + parent->genome[i]) / 2.0;
            else if (crossover_method == CROSSOVER_SUBSTITUTION)
                child->genome[i] = parent->genome[i];
        }
    }
    else
    {
        // copy left side of crossover point
        for (int i = pos; i >= 0; i--)
        {
            if (crossover_method == CROSSOVER_AVERAGE)
                child->genome[i] = (child->genome[i] + parent->genome[i]) / 2.0;
            else if (crossover_method == CROSSOVER_SUBSTITUTION)
                child->genome[i] = parent->genome[i];
        }
    }
}

// performs crossover step on supplied population
void perform_crossover(chromozome** population, int popsize)
{
    // suitable parents
    std::vector<chromozome*> parent_list;

    // get lowest and highest fitness
    double fitmax = population[0]->fitness;
    double fitmin = population[0]->fitness;
    for (int i = 0; i < population_size; i++)
    {
        if (population[i]->fitness > fitmax)
            fitmax = population[i]->fitness;
        if (population[i]->fitness < fitmin)
            fitmin = population[i]->fitness;
    }

    // get the fitness distance
    double distance = fitmax - fitmin;
    double chance;

    int pos = 0;

    // chance to crossover - BASE*((fitness - fitmin)/distance)
    for (int i = 0; i < popsize; i++)
    {
        // lower crossover rate for chromosomes with lower fitness
        chance = crossover_base_chance_pct * ((population[i]->fitness - fitmin) / distance);
        // roll the dice
        if (frand() < chance)
            parent_list.push_back(population[i]);
    }

    // if no suitable parents were selected, there's nothing to do
    if (parent_list.size() == 0)
        return;

    // prepare parent array from vector
    chromozome** parent_array = new chromozome*[parent_list.size()];
    for (int i = 0; i < (int) parent_list.size(); i++)
        parent_array[i] = parent_list[i];

    chromozome* tmpsel;

    if (parent_list.size() == 1)
    {
        // 1 element list, cross the only one chromosome with random from whole population
        do
        {
            tmpsel = select_random_chromozome(population, popsize, parent_list[0]);
        } while (tmpsel == parent_list[0] && tmpsel != NULL);

        // if some chromosome was selected...
        if (tmpsel != NULL)
        {
            // .. decide which parent has greated fitness, and consider the one with lower as "child"
            // therefore rewrite his genome
            if (parent_list[0]->fitness > tmpsel->fitness)
                crossover_single(parent_list[0], tmpsel);
            else
                crossover_single(tmpsel, parent_list[0]);
        }
    }
    else
    {
        // more elements list, cross between parent list members
        for (int i = 0; i < (int) parent_list.size(); i++)
        {
            do
            {
                tmpsel = select_random_chromozome(parent_array, parent_list.size(), parent_list[i]);
            } while (tmpsel == parent_list[i] && tmpsel != NULL);

            // if some chromosome was selected...
            if (tmpsel != NULL)
            {
                // .. decide which parent has greated fitness, and consider the one with lower as "child"
                // therefore rewrite his genome
                if (parent_list[i]->fitness > tmpsel->fitness)
                    crossover_single(parent_list[i], tmpsel);
                else
                    crossover_single(tmpsel, parent_list[i]);
            }
        }
    }

    // clean things up!
    parent_list.clear();
    delete[] parent_array;
}

// performs mutation within one chromosome
void mutation_single(chromozome* chrom)
{
    chrom->genome[rand() % genome_size] = parameter_limit_min + (frand() * (double)(parameter_limit_max - parameter_limit_min));
}

// proceeds mutation on supplied population
void perform_mutation(chromozome** population)
{
    // go through all chromosomes in population
    for (int i = 0; i < population_size; i++)
    {
        // and roll the dice for the mutation
        if (frand() < mutation_base_chance_pct)
            mutation_single(population[i]);
    }
}

// bubble sort...
void pop_sort(chromozome** population, int popsize)
{
    chromozome* tmp;

    for (int i = 0; i < popsize; i++)
    {
        for (int j = 0; j < popsize - 1; j++)
        {
            if (population[j]->fitness > population[j + 1]->fitness)
            {
                tmp = population[j];
                population[j] = population[j + 1];
                population[j + 1] = tmp;
            }
        }
    }
}

// proceeds genocide on supplied population
void perform_genocide(chromozome** population)
{
    // sort from least fit to most fit
    pop_sort(population, population_size);

    // "kill" one third of population
    for (int i = 0; i < (population_size / 3); i++)
    {
        // renew their genome
        for (int j = 0; j < genome_size; j++)
            population[i]->genome[j] = parameter_limit_min + (frand() * (double)(parameter_limit_max - parameter_limit_min));
    }
}

// prints out genome of supplied chromosome with its fitness
void print_genome(chromozome* chrom)
{
    cout << "Genome: [";
    for (int j = 0; j < genome_size; j++)
    {
        cout << chrom->genome[j];
        if (j == genome_size - 1)
            cout << "]";

        cout << ", ";
    }
    cout << "fitness: " << chrom->fitness << endl;
}

// prints out information about population
void print_population(int generation, chromozome** population)
{
    // calculate average fitness
    double sum = 0.0;
    for (int i = 0; i < population_size; i++)
        sum += population[i]->fitness;

    cout << "Generation: " << generation << ", average fitness: " << (sum/(double)population_size) << endl;
}

// retrieves the fittest one chromosome from supplied population
chromozome* get_fittest(chromozome** population)
{
    chromozome* fittest = NULL;

    for (int i = 0; i < population_size; i++)
    {
        if (fittest == NULL || fittest->fitness < population[i]->fitness)
            fittest = population[i];
    }

    return fittest;
}

// returns the fittest chromosome back into population, substitutes the least fit one
void return_fittest(chromozome** population, chromozome* fittest)
{
    int least = 0;

    for (int i = 0; i < population_size; i++)
    {
        if (population[least]->fitness > population[i]->fitness)
            least = i;
    }

    double* tmpptr = population[least]->genome;

    memcpy(population[least], fittest, sizeof(chromozome));
    population[least]->genome = tmpptr;
    memcpy(population[least]->genome, fittest->genome, sizeof(double)*genome_size);
}

// is variable already present in variable vector?
bool is_var_present(char var)
{
    for (std::list<char>::iterator itr = variable_map.begin(); itr != variable_map.end(); ++itr)
    {
        if (*itr == var)
            return true;
    }

    return false;
}

int main(int argc, char** argv)
{
    // 7 parameters are mandatory
    if (argc < 7)
    {
        cout << "Parameters: " << endl
            << " <left_side_expression>" << endl
            << " <right_side_constant>" << endl
            << " <population_size>" << endl
            << " <generation_count>" << endl
            << " <variable_limit_down>" << endl
            << " <variable_limit_up>" << endl
            << " [crossover_chance = 0.3]" << endl
            << " [mutation_chance = 0.1]" << endl
            << " [crossover_method = A]" << endl
            << " [genocides = 0]" << endl;

        cout << endl << "Crossover method: A = average, S = substitution" << endl;
        return 1;
    }

    // load parameters
    char* expression = argv[1];
    rightside = atof(argv[2]);
    population_size = atoi(argv[3]);
    generation_count_limit = atoi(argv[4]);
    parameter_limit_min = atof(argv[5]);
    parameter_limit_max = atof(argv[6]);
    if (argc >= 8)
        crossover_base_chance_pct = atof(argv[7]);
    if (argc >= 9)
        mutation_base_chance_pct = atof(argv[8]);
    if (argc >= 10)
    {
        char method = argv[9][0];
        if (method == 'A')
            crossover_method = CROSSOVER_AVERAGE;
        else if (method == 'S')
            crossover_method = CROSSOVER_SUBSTITUTION;
        else
        {
            cout << "Unsupported crossover method, falling back to averaging" << endl;
            crossover_method = CROSSOVER_AVERAGE;
        }
    }
    int genocides = 0;
    if (argc >= 11)
        genocides = atoi(argv[10]);

    int err;
    char* err_ptr;

    // parse input expression
    exp_stack = sy_generate_rpn_stack(argv[1], &err, &err_ptr);
    if (exp_stack == NULL)
    {
        cout << "You have an error in supplied expression!" << endl;
        return 1;
    }

    // prepare evaluation
    int varcount = 0;
    int hvar = 0;
    rpn_element* el;
    for (int i = 0; i < exp_stack->size; i++)
    {
        el = (rpn_element*)stck_get(exp_stack, i);
        if (!el)
            continue;
        // add all variables (just once each) to variable list
        if (el->type == RPN_TOKEN_VARIABLE)
        {
            // every variable just once
            if (!is_var_present(el->value.as_variable))
            {
                variable_map.push_back(el->value.as_variable);
                if (el->value.as_variable > hvar)
                    hvar = el->value.as_variable;
                varcount++;
            }
        }
    }

    // genome size is variable count
    genome_size = varcount;
    // init evaluation map
    eval_map = new double[hvar+1];

    // init PRNG
    srand((unsigned int) time(NULL));

    // init algorithm

    int generation = 0;
    // allocate what we need
    chromozome** population = new chromozome*[population_size];
    chromozome* fittest = new chromozome();
    chromozome* tmpfit;

    // generate random first population
    for (int i = 0; i < population_size; i++)
    {
        population[i] = new chromozome();
        population[i]->genome = new double[genome_size];

        for (int j = 0; j < genome_size; j++)
            population[i]->genome[j] = parameter_limit_min + (frand() * (double)(parameter_limit_max - parameter_limit_min));

        population[i]->fitness = fitness_function(population[i]);
    }

    // store the fittest one for later reuse (elitism strategy)
    tmpfit = get_fittest(population);
    memcpy(fittest, tmpfit, sizeof(chromozome));
    fittest->genome = new double[genome_size];
    memcpy(fittest->genome, tmpfit->genome, sizeof(double)*genome_size);

    // while we can evaluate next generation..
    while (++generation <= generation_count_limit)
    {
        // if it's time to genocide (if set, the genocides are distributed uniformly in generations)
        if (genocides > 0 && generation == (int)(generation_count_limit / (genocides+1)))
        {
            perform_genocide(population);
            recalculate_fitness(population, population_size);
            // the fittest one always survives
            return_fittest(population, fittest);
        }

        // at first, perform crossover
        perform_crossover(population, population_size);
        // then mutation
        perform_mutation(population);

        // return the old fittest on back - he will kick out the least fit one
        return_fittest(population, fittest);

        // get the new fittest one
        tmpfit = get_fittest(population);
        if (tmpfit->fitness > fittest->fitness)
        {
            double* tmpptr = fittest->genome;
            memcpy(fittest, tmpfit, sizeof(chromozome));
            fittest->genome = tmpptr;
            memcpy(fittest->genome, tmpfit->genome, sizeof(double)*genome_size);
        }

        // update fitness of all chromosomes
        recalculate_fitness(population, population_size);
        // and print something nice about this generation
        print_population(generation, population);
    }

    cout << "----------" << endl;
    cout << "Expression: " << argv[1] << " = " << rightside << endl;

    cout << "Fittest: " << endl;
    print_genome(fittest);
    cout << endl << "Deviation: " << 1.0/fittest->fitness << endl;

    // clean up all the mess we've done
    delete[] eval_map;
    for (int i = 0; i < population_size; i++)
    {
        delete[] population[i]->genome;
        delete population[i];
    }
    delete[] population;
    delete[] fittest->genome;
    delete fittest;

    stck_destroy(exp_stack);

    return 0;
}

