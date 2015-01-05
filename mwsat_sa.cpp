#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <algorithm>

struct inst
{
	std::vector<int> weights;
	std::vector<std::vector<int> > formula;
	int id;
	int num_clauses;//length of formula outer vector
	int num_vars;//length of weights
};
struct soln
{
	uint64_t config;
	int id;
	int num_items;
	int cost;
};
struct SAParams
{
	double initTemp;
	double finalTemp;
	int numIterations;
	double geoCoolFactor;
	bool saveProgress;
};

std::ofstream progressFile;

//helpers
void printInst(inst instance)
{
	printf("num_clauses: %d, num_vars: %d\n", instance.num_clauses, instance.num_vars);
	printf("Formula:\n");
	for (int i = 0; i < instance.num_clauses; i++)
	{
		printf("(");
		for (int j = 0; j < instance.formula[i].size(); ++j)
		{
			printf("x%d", std::abs(instance.formula[i][j]));
			if (instance.formula[i][j] < 0) printf("'");
			if (j < instance.formula[i].size() - 1) printf(" + ");
		}
		if (i < instance.num_clauses - 1) printf(").");
		else printf(")");
	}
	printf("\n");
}

void printSoln(soln solution)
{
	printf("Id: %d, n: %d, cost: %d, config: %lld\n", solution.id, solution.num_items, solution.cost, solution.config);
}

int GetInstances(std::vector<inst> &instances, const char* filename, bool verbose) 
{
	std::string line;
  	std::ifstream infile(filename);
  	std::string tempstr;
  	int tempint;

  	if (infile)
  	{
	  	while (std::getline(infile, line))
	  	{
  			inst newInst;
  			std::istringstream iss(line);

	  		iss >> tempstr;
	  		if (tempstr.compare(std::string("cnf")) == 0)
	  		{
	  			iss >> newInst.num_vars;
	  			iss >> newInst.num_clauses;

	  			std::getline(infile, line);
	  			std::istringstream iss(line);

	  			iss >> tempstr;
	  			if (tempstr.compare(std::string("w")) != 0) return 1;

	  			for (int i = 0; i < newInst.num_vars; ++i)
	  			{
	  				iss >> tempint;
	  				newInst.weights.push_back(tempint);
	  			}

	  			for (int i = 0; i < newInst.num_clauses; ++i)
	  			{
	  				newInst.formula.push_back(std::vector<int>());

	  				std::getline(infile, line);
	  				std::istringstream iss(line);

	  				int clauseSize;
	  				iss >> clauseSize;

	  				for (int j = 0; j < clauseSize; ++j)
	  				{
	  					iss >> tempint;
	  					newInst.formula[i].push_back(tempint);
	  				}
	  			}
	  		}
	  		else
	  		{
	  			return 1;
	  		}

	  		instances.push_back(newInst);

	  		if (verbose)
	  		{
	  			printf("Pushed New Instance:\n");
	  			printInst(newInst);
	  		}
	  	}
  	}

  	return 0;
}

std::vector<bool> GenerateInitState(bool isRandom, int n)
{
	std::vector<bool> state;
	if(isRandom)
	{
		for (int i = 0; i < n; ++i)
		{
			state.push_back((bool)(rand() % 2));
		}
	}
	else
	{
		//TODO: greedy to start with local minimum
	}
	return state;
}

std::vector<bool> RandomNeighbour(std::vector<bool> state)
{
	int randIndex = rand() % state.size();
	state[randIndex] = !state[randIndex];

	return state;
}

int Weight(std::vector<bool> state, inst instance)
{
	int stateWeight = 0;
	for (int i = 0; i < state.size(); ++i)
	{
		stateWeight += instance.weights[i] * (int)state[i];
	}
	return stateWeight;
}

bool Feasible(std::vector<bool> state, inst instance)//check if satisfied
{
	bool outer = 1;//start with 1 since it will be AND'd with clauses (any 0 -> formula unsatisfied)
	for (int i = 0; i < instance.num_clauses; ++i)
	{
		bool inner = 0;//start with 0 since it will be OR'd with variables (any 1 -> clause satisfied)
		for (int j = 0; j < instance.formula[i].size(); ++j)
		{
			if (instance.formula[i][j] < 0)
			{
				inner = inner | !(state[std::abs(instance.formula[i][j]) - 1]);
			}
			else
			{
				inner = inner | state[instance.formula[i][j] - 1];
			}
		}
		outer = outer & inner;
	}

	return outer;
}

void Repair(std::vector<bool> &state, inst instance)
{
	while(!Feasible(state, instance))
	{
		//remove random item
		state = RandomNeighbour(state);
	}
}

bool Accept(int weightDiff, int t)
{
	if (t == 0) return 0;
	bool isAccept = (double)rand() / RAND_MAX < exp((weightDiff)/t);
	return isAccept;
}

//Simulated Annealing Algorithm
int SimAnneal(inst instance, SAParams params)
{
	int iterCount = 0;
	int highestWeight = 0;

	std::clock_t start;
    double duration;
    start = std::clock();

	double t = params.initTemp;
	std::vector<bool> state = GenerateInitState(true, instance.num_vars);

	while(t > params.finalTemp)
	{
		for (int i = 0; i < params.numIterations; ++i)
		{
			iterCount++;
			std::vector<bool> newState = RandomNeighbour(state);
			Repair(newState, instance);//checks feasibility and repairs if necessary.
			int curWeight = Weight(newState, instance);
			highestWeight = curWeight > highestWeight ? curWeight : highestWeight;//save highest cost
			int weightDiff = curWeight - Weight(state, instance);
			if (weightDiff > 0)
			{
				state = newState;
			}
			else if (Accept(weightDiff, t))
			{
				state = newState;
			}

			if (params.saveProgress)
			{
				progressFile << iterCount << ", " << Weight(state, instance) << "\n";
			}
		}
		t = t * params.geoCoolFactor;
	}

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return highestWeight;
}

int SimAnnealWrapper(std::vector<inst> instances, std::vector<soln> &solutions, SAParams params)
{
	while (!instances.empty())
	{
		inst curInst = instances.back();
		instances.pop_back();

		if (params.saveProgress)
		{
			progressFile.open("progress.csv");
		}
		int cost = SimAnneal(curInst, params);
		printf("%d\n", cost);
		if (params.saveProgress)
		{
			progressFile.close();
		}
	}

	return 0;
}


int main(int argc, const char** argv) 
{
	if (argc < 7 || argc > 8){
		printf("Error: Unexpected argument count.\n");
		printf("Usage: <filename> <initial temp> <final temp> <num iterations> <cooling factor> <output progress?> [-v]\n\n"); 
		return 1;
	}

	SAParams params = {
		.initTemp = atof(argv[2]),
		.finalTemp = atof(argv[3]),
		.numIterations = atoi(argv[4]),
		.geoCoolFactor = atof(argv[5]),
		.saveProgress = (bool)atoi(argv[6])
	};

	if (params.saveProgress)
	{
		printf("Saving Progress...\n");
	}

	bool verbose;
	if (argc == 8)
	{
		verbose = (std::strncmp(argv[7], "-v", 2) == 0);
	}

	srand(time(NULL));//seed the random number generator

	std::vector<inst> instances;
	std::vector<soln> solutions;

	GetInstances(instances, argv[1], verbose);

	SimAnnealWrapper(instances, solutions, params);

	instances.clear();
	solutions.clear();

	return 0;
}



