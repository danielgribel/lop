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
const double epslon = 0.0001;

//int** simDegrees;
//int* pos;

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

	//pos[v[i]] = j;
	//pos[v[j]] = i;
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
		//pos[v[k+1]] = k;
	}
	v[j] = tmp;
	//pos[tmp] = i;
}

// for j < i
void insert_(int *v, const int i, const int j) {
	int tmp = v[i];
	for (int k = i; k > j; k--) {
		v[k] = v[k-1];
		//pos[v[k-1]] = k;
	}
	v[j] = tmp;
	//pos[tmp] = i;
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
    	if(randItems[i] < randItems[i+1]) {
    		cost = testSwap(s, randItems[i], randItems[i+1]);	
    		swap(s, randItems[i], randItems[i+1]);
    	} else {
    		cost = testSwap(s, randItems[i+1], randItems[i]);	
    		swap(s, randItems[i+1], randItems[i]);
    	}
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

int testInsertConsec(int* v, int i, int j) {
	return diff[v[j]][v[i]];
}

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
				//delta = delta + testInsertConsec(s, i1, j);
				delta = delta + diff[s[i1]][s[j]];
				if(delta > 0) {
					insert_(s, i1, j);
					improving = true;
					delta = 0;
					i1 = j;
				}
			}
			delta = 0;
			for(int j = i+1; j < n; j++) {
				//delta = delta + testInsertConsec(s, j, i1);
				delta = delta + diff[s[j]][s[i1]];
				if(delta > 0) {
					insert(s, i1, j);
					improving = true;
					delta = 0;
					i1 = j;
				}
			}
			//checkSwap(s, i);
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

    //pos = new int[n];
    
    for(int i = k-1; i >= 0; --i) {
        ki = q.top().second;
        n_closest[i] = ki;
        q.pop();
        //pos[ki] = i;
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
	
	/*simDegrees = new int*[n];
	const int degreesSize = 4;

	for(int i = 0; i < n; i++) {
		simDegrees[i] = new int[degreesSize];
	}

	for(int i = 0; i < n; i++) {
		if(i == 0) {
			simDegrees[kSmallest[i]][0] = kSmallest[i+1];
			simDegrees[kSmallest[i]][1] = kSmallest[i+2];
			simDegrees[kSmallest[i]][2] = kSmallest[i+3];
			simDegrees[kSmallest[i]][3] = kSmallest[i+4];
		} else if(i == 1) {
			simDegrees[kSmallest[i]][0] = kSmallest[i-1];
			simDegrees[kSmallest[i]][1] = kSmallest[i+1];
			simDegrees[kSmallest[i]][2] = kSmallest[i+2];
			simDegrees[kSmallest[i]][3] = kSmallest[i+3];
		} else if(i == n-1) {
			simDegrees[kSmallest[i]][0] = kSmallest[i-1];
			simDegrees[kSmallest[i]][1] = kSmallest[i-2];
			simDegrees[kSmallest[i]][2] = kSmallest[i-3];
			simDegrees[kSmallest[i]][3] = kSmallest[i-4];
		} else if(i == n-2) {
			simDegrees[kSmallest[i]][0] = kSmallest[i+1];
			simDegrees[kSmallest[i]][1] = kSmallest[i-1];
			simDegrees[kSmallest[i]][2] = kSmallest[i-2];
			simDegrees[kSmallest[i]][3] = kSmallest[i-3];
		} else {
			simDegrees[kSmallest[i]][0] = kSmallest[i-1];
			simDegrees[kSmallest[i]][1] = kSmallest[i-2];
			simDegrees[kSmallest[i]][2] = kSmallest[i+1];
			simDegrees[kSmallest[i]][3] = kSmallest[i+2];
		}
	}*/

	delete [] sumOut;
	delete [] sumIn;
	delete [] delta;

	return kSmallest;
}

void verifySolution(std::string inputFile, int* bestSolution, int bestCost, double elapsedSecs) {
	int* sol = new int [n];
	for(int i = 0; i < n; i++) {
		sol[i] = bestSolution[i];
	}
	int c = getCost(sol, matrix, n);
	avgCost = avgCost + c;

	std::cout << inputFile << " ";
	printf("%d ", c);
	printf("%d ", c == bestCost);
	printf("%f\n", elapsedSecs);
	
	delete [] sol;
}

void run(std::string inputFile) {
	srand(1607);
	loadData(inputFile);

	clock_t begin = clock();

	int bestCost = 0;
	int* bestSolution;
	int it = 1;
	int lastImprovement = 0;
	const int numExchanges = 10;
	const int itNoImprovement = 100;
	const int maxIterations = 500;
	
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
		if((1.0*cost) < ((1-epslon)*bestCost)) {
			s = bestSolution;
		}
		//printf("f(%d) = %d lastImprovement = %d\n", it, bestCost, lastImprovement);
		it++;
	}

	double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
	
	avgTime = avgTime + elapsedSecs;

	verifySolution(inputFile, bestSolution, bestCost, elapsedSecs);

	deleteMatrix(matrix, n);
	deleteMatrix(diff, n);
	//deleteMatrix(simDegrees, n);
	//delete [] pos;
	delete [] s;	
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
	const char* FOLDER = "data/n2000/";

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