#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <ctime>
#include <cstdlib>
#include <set>
#include <queue>
#include <functional>
#include <limits>

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();
const string DATASET = "data/n0500/n0500d100-5";
int** matrix;
int cost;
int n;

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

int testInsert(int *v, const int i, const int j) {
	int newcost = cost;

	for(int k = i+1; k <= j; k++) {
		newcost = newcost - matrix[v[i]][v[k]] + matrix[v[k]][v[i]];
	}

	return newcost;
}

void insert(int *v, const int i, const int j) {
	int tmp = v[i];
	for (int k = i; k < j; k++) {
		v[k] = v[k+1];
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
    int randItems[n];
    int a[k];
    int a_[k];

    for(int i = 0; i < n; i++) {
        s[i] = s0[i];
    }

    shuffle(randItems, n);
    shuffle(a_, k);

    for(int i = 0; i < k; i++) {
      a[i] = randItems[i];
    }

    int temp[k];

    for(int i = 0; i < k; i++) {
      temp[i] = s[a[a_[i]]];
    }

    for(int i = 0; i < k; i++) {
      s[a[i]] = temp[i];
    }

    return s;
}

void loadData(string input) {
	int x, y;
	ifstream in(input.c_str());

	if(!in) {
		cout << "Cannot open file.\n";
		return;
	}

	in >> n;
	
	matrix = new int*[n];
  
	for(int i = 0; i < n; i++) {
		matrix[i] = new int[n];
	}

	for(x = 0; x < n; x++) {
		for(y = 0; y < n; y++) {
			in >> matrix[x][y];
		}
	}

	in.close();
}

void localSearch(int* s, const int n) {
	//int sh[n];
	//shuffle(sh, n);

	int i;
	int it = 0;
	int newcost;
	bool improving = true;
	
	while(improving) {
		it++;
		for(int i = 0; i < n; i++) {
			//i = sh[i1];
			improving = false;
			/*for(int j = i+1; j < n; j++) {
				newcost = testInsert(s, i, j);
				if(newcost > cost) {
					insert(s, i, j);
					cost = newcost;
					improving = true;
				}
			}*/
			newcost = testRotateLeft(s, i, n);
			if(newcost > cost) {
				rotateLeft(s, i, n);
				cost = newcost;
				improving = true;
			}
			//if(!improving) {
			for(int j = i+1; j < n; j++) {
				newcost = testSwap(s, i, j);
				if(newcost > cost) {
					swap(s, i, j);
					cost = newcost;
					improving = true;
				}
			}
			//}
		}
	}
}

int main() {
	clock_t begin = clock();

	srand(1607);
	loadData(DATASET);

	int s0[n];
	shuffle(s0, n);

	int bestCost = 0;
	int* bestSolution = new int[n];
	int p = 1;
	int numPert = 100;
	int numExchanges = 5;

	cost = getCost(s0, matrix, n);

	localSearch(s0, n);
	printf("f(0) = %d\n", cost);

	if(cost > bestCost) {
		bestCost = cost;
		bestSolution = s0;
	}

	int* s = s0;

	while(p <= numPert) {
		s = perturbation(s, numExchanges, n);
		cost = getCost(s, matrix, n);
		localSearch(s, n);		
		if(cost > bestCost) {
			bestCost = cost;
			bestSolution = s;
		}
		printf("f(%d) = %d\n", p, bestCost);
		p++;
	}

	//printf("g(%d) = %.15g\n", i, bestCost);
	
	int sol[n];

	for(int i = 0; i < n; i++) {
		s[i] = bestSolution[i];
	}

	int c = getCost(s, matrix, n);
	printf("c = %d\n", c);

	double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
	cout << elapsedSecs << endl;

	return 0;
}