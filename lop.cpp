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

#define pii pair<int, int>
#define pdi pair<double, int>

const double MAX_FLOAT = std::numeric_limits<double>::max();
const string DATASET = "data/n1000/n1000d050-4";
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

// for i < j
int testSwapConsec(int *v, const int i, const int j) {
	int delta = matrix[v[j]][v[i]] - matrix[v[i]][v[j]];
	return delta;
}

// for j < i
int testSwapConsec_(int *v, const int i, const int j) {
	int delta = matrix[v[i]][v[j]] - matrix[v[j]][v[i]];
	return delta;
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
			for(int j = i-1; j >= 0; j--) {
				newcost = testInsert_(s, i, j);
				if(newcost > cost) {
					insert_(s, i, j);
					cost = newcost;
					improving = true;
				}
			}
			for(int j = i+1; j < n; j++) {
				newcost = testInsert(s, i, j);
				if(newcost > cost) {
					insert(s, i, j);
					cost = newcost;
					improving = true;
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
		for(int i = 0; i < n; i++) {
			//i = sh[i0];
			delta = 0;
			i1 = i;
			for(int j = i-1; j >= 0; j--) {
				delta = delta + testSwapConsec_(s, i1, j);
				if(delta > 0) {
					insert_(s, i1, j);
					improving = true;
					delta = 0;
					i1 = j;
				}
			}
			delta = 0;
			for(int j = i+1; j < n; j++) {
				delta = delta + testSwapConsec(s, i1, j);
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
    priority_queue< pii, vector < pii >, greater< pii > > q;

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

int* degrees() {
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
		}
		delta[i] = sumOut[i] - sumIn[i];
	}

	int* kSmallest = kSmallestIndices(delta, n, n);

	delete [] sumOut;
	delete [] sumIn;
	delete [] delta;
	
	return kSmallest;
}

void verifySolution(int* bestSolution, int bestCost) {
	int* sol = new int [n];
	for(int i = 0; i < n; i++) {
		sol[i] = bestSolution[i];
	}
	int c = getCost(sol, matrix, n);
	printf("cost = %d\n", c);
	if(c == bestCost) {
		printf("Assert = %c\n", 'T');	
	} else {
		printf("Assert = %c\n", 'F');
	}
	delete [] sol;
}

int main() {
	clock_t begin = clock();

	srand(1607);
	
	int* s = new int[n];
	int bestCost = 0;
	int* bestSolution;
	int p = 1;
	const int numPert = 100;
	const int numExchanges = 5;

	loadData(DATASET);
	s = degrees();
	cost = getCost(s, matrix, n);
	localSearch(s, n);
	printf("f(0) = %d\n", cost);
	bestCost = cost;
	bestSolution = s;

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

	verifySolution(bestSolution, bestCost);

	double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
	cout << elapsedSecs << endl;

	delete [] s;

	return 0;
}