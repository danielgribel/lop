#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <queue>
#include <functional>
#include <limits>
#include <dirent.h>
#include <vector>
#include <algorithm>

#define pii std::pair<int, int>
#define pdi std::pair<double, int>

const double MAX_FLOAT = std::numeric_limits<double>::max();
int** matrix;
int** diff;
int cost;
int n;
double avgCost;
double avgTime;
const double epslon = 0.0000000;

int getCost(int* s, int** matrix, int n) {
	int cost = 0;
	int i1;
	int j1;
	for(int i = 0; i < n-1; i++) {
		i1 = s[i];
		for(int j = i+1; j < n; j++) {
			j1 = s[j];
			cost = cost + matrix[i1][j1];
		}
	}
	return cost;
}

void print(const int *v, const int size) {
	if (v != 0) {
		for (int i = 0; i < size; i++) {
			printf("%d ", v[i]);
		}
		printf("\n");
	}
}

int testSwap(int *v, const int i, const int j) {
	int firstElement = v[i];
	int secondElement = v[j];
	int newcost = cost;

	newcost = newcost - matrix[firstElement][secondElement];
	newcost = newcost + matrix[secondElement][firstElement];

	for(int i1 = i+1; i1 < j; i1++) {
		newcost = newcost - matrix[firstElement][v[i1]] + matrix[secondElement][v[i1]];
		newcost = newcost - matrix[v[i1]][secondElement] + matrix[v[i1]][firstElement];
	}

	return newcost;
}

void swap(int *v, const int i, const int j) {
	int t;
	t = v[i];
	v[i] = v[j];
	v[j] = t;
}

int testRotateLeft(int *v, const int i, const int n) {
	int newcost = cost;

	for(int j = i+1; j < n; j++) {
		newcost = newcost - matrix[v[i]][v[j]] + matrix[v[j]][v[i]];
	}

	return newcost;
}

void rotateLeft(int *v, const int start, const int n) {
	int tmp = v[start];
	for (int i = start; i < n-1; i++) {
		v[i] = v[i+1];
	}
	v[n-1] = tmp;
}

// for i < j
int testInsert(int *v, const int i, const int j) {
	int newcost = cost;

	for(int k = i+1; k <= j; k++) {
		newcost = newcost - matrix[v[i]][v[k]] + matrix[v[k]][v[i]];
	}

	return newcost;
}

// for j < i
int testInsert_(int *v, const int i, const int j) {
	int newcost = cost;

	for(int k = j; k <= i; k++) {
		newcost = newcost - matrix[v[k]][v[i]] + matrix[v[i]][v[k]];
	}

	return newcost;
}

// for i < j
void insert(int *v, const int i, const int j) {
	int tmp = v[i];
	for (int k = i; k < j; k++) {
		v[k] = v[k+1];
	}
	v[j] = tmp;
}

// for j < i
void insert_(int *v, const int i, const int j) {
	int tmp = v[i];
	for (int k = i; k > j; k--) {
		v[k] = v[k-1];
	}
	v[j] = tmp;
}

void init(int *v, int n) {
	for(int i = 0; i < n; i++) {
		v[i] = i;
	}
}

void shuffle(int *myArray, size_t n) {
    for(int i = 0; i < n; i++) {
        myArray[i] = i;
    }
    if(n > 1) {
        size_t i;
        for (int i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = myArray[j];
            myArray[j] = myArray[i];
            myArray[i] = t;
        }
    }
}

int* perturbation(int* s0, int k, int n) {
	int* s = new int[n];
    int* randItems = new int[n];

    shuffle(randItems, n);

    for(int i = 0; i < n; i++) {
    	s[i] = s0[i];
    }

    for(int i = 0; i < 2*k; i = i+2) {
    	cost = testSwap(s, randItems[i], randItems[i+1]);
        swap(s, randItems[i], randItems[i+1]);
    }

    delete [] randItems;

    return s;
}

void loadData(std::string input) {
	int x, y;
	std::ifstream in(input.c_str());

	if(!in) {
		std::cout << "Cannot open file.\n";
		return;
	}

	in >> n;
	
	matrix = new int*[n];
  	diff = new int*[n];

	for(int i = 0; i < n; i++) {
		matrix[i] = new int[n];
		diff[i] = new int[n];
	}

	for(x = 0; x < n; x++) {
		for(y = 0; y < n; y++) {
			in >> matrix[x][y];
		}
	}

	in.close();
}

void deleteMatrix(int** matrix, int n) {
	for(int i = 0; i < n; i++) {
		delete [] matrix[i];
	}
	delete [] matrix;
}

/*void localSearch(int* s, const int n) {
	int* sh = new int[n];
	shuffle(sh, n);

	int i, j;
	int it = 0;
	int newcost;
	bool improving = true;
	
	while(improving) {
		it++;
		for(int i = 0; i < n; i++) {
			improving = false;
			for(int j = i+1; j < n; j++) {
				newcost = testInsert(s, i, j);
				if(newcost > cost) {
					insert(s, i, j);
					cost = newcost;
					improving = true;
				}
			}
			newcost = testRotateLeft(s, i, n);
			if(newcost > cost) {
				rotateLeft(s, i, n);
				cost = newcost;
				improving = true;
			}
			if(!improving) {
				for(int j = i+1; j < n; j++) {
					newcost = testSwap(s, i, j);
					if(newcost > cost) {
						swap(s, i, j);
						cost = newcost;
						improving = true;
					}
				}
			}
		}
	}

	delete [] sh;
}*/

void localSearch(int* s, const int n) {
	int* sh = new int[n];
	shuffle(sh, n);

	int i, j, i1;
	int it = 0;
	int delta;
	bool improving = true;

	while(improving) {
		it++;
		improving = false;
		for(int i0 = 0; i0 < n; i0++) {
			i = sh[i0];
			delta = 0;
			i1 = i;
			for(int j = i-1; j >= 0; j--) {
				delta = delta + diff[s[i1]][s[j]];
				//if((1.0*(cost+delta)) > ((1-epslon)*cost)) {
				if(delta > 0) {
					insert_(s, i1, j);
					improving = true;
					delta = 0;
					i1 = j;
				}
			}
			delta = 0;
			for(int j = i+1; j < n; j++) {
				delta = delta + diff[s[j]][s[i1]];
				//if((1.0*(cost+delta)) > ((1-epslon)*cost)) {
				if(delta > 0) {
					insert(s, i1, j);
					improving = true;
					delta = 0;
					i1 = j;
				}
			}
		}
	}
	cost = getCost(s, matrix, n);
	delete [] sh;
}

int* kSmallestIndices(int* v, int n, int k) {
    int* n_closest = new int[k];
    std::priority_queue< pii, std::vector < pii >, std::greater< pii > > q;

    for(int i = 0; i < n; ++i) {
        q.push(pii(v[i], i));
    }
    
    int ki;
    
    for(int i = k-1; i >= 0; --i) {
        ki = q.top().second;
        n_closest[i] = ki;
        q.pop();
    }

    return n_closest;
}

// initial solution based on the difference between out-degree and in-degree
int* initialSolution() {
	int* sumOut = new int[n];
	int* sumIn = new int[n];
	int* delta = new int[n];
	
	for(int i = 0; i < n; i++) {
		sumOut[i] = 0;
		sumIn[i] = 0;
	}

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			sumOut[i] = sumOut[i] + matrix[i][j];
			sumIn[i] = sumIn[i] + matrix[j][i];
			diff[i][j] = matrix[i][j] - matrix[j][i];
		}
		delta[i] = sumOut[i] - sumIn[i];
	}

	int* kSmallest = kSmallestIndices(delta, n, n);

	delete [] sumOut;
	delete [] sumIn;
	delete [] delta;
	
	return kSmallest;
}

void verifySolution(std::string inputFile, int* bestSolution, int bestCost) {
	int* sol = new int [n];
	for(int i = 0; i < n; i++) {
		sol[i] = bestSolution[i];
	}
	int c = getCost(sol, matrix, n);
	avgCost = avgCost + c;

	std::cout << inputFile << " ";
	printf("%10d ", c);
	printf("%5d ", c == bestCost);
	
	delete [] sol;
}

void run(std::string inputFile) {
	clock_t begin = clock();
	
	srand(1607);

	int bestCost = 0;
	int* bestSolution;
	int it = 1;
	int lastImprovement = 0;
	const int numExchanges = 5;
	const int itNoImprovement = 100;
	const int maxIterations = 500;
	
	loadData(inputFile);
	int* s = initialSolution();
	localSearch(s, n);
	//printf("f(0) = %d\n", cost);
	bestCost = cost;
	bestSolution = s;

	while((it-lastImprovement < itNoImprovement) && (it < maxIterations)) {
		s = perturbation(s, numExchanges, n);
		localSearch(s, n);
		if(cost > bestCost) {
			bestCost = cost;
			bestSolution = s;
			lastImprovement = it;
		}
		//printf("f(%d) = %d lastImprovement = %d\n", it, bestCost, lastImprovement);
		it++;
	}

	verifySolution(inputFile, bestSolution, bestCost);

	double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
	printf("%10f\n", elapsedSecs);

	avgTime = avgTime + elapsedSecs;

	delete [] s;
	deleteMatrix(matrix, n);
	deleteMatrix(diff, n);
}

std::vector<std::string> listFiles(const char* folder) {
    DIR *d;
    struct dirent *dir;
    int i = 0;
    std::vector<std::string> files;

    d = opendir(folder);

    if(d) {
        while((dir = readdir(d)) != NULL) {
            i++;
            if(dir->d_name[0] != '.') {
                files.push_back(dir->d_name);    
            }
        }
        closedir(d);
    }
    std::sort( files.begin(), files.end() );
    return files;
}

int main() {
	const char* FOLDER = "data/n0500/";

	avgCost = 0.0;
	avgTime = 0.0;

    std::vector<std::string> files = listFiles(FOLDER);

    for(int i = 0; i < files.size(); i++) {
    	run(FOLDER + files[i]);
    }
    
    printf("avg cost = %.10f\n", avgCost/files.size());
    printf("avg time = %.10f\n", avgTime/files.size());
	
	return 0;
}