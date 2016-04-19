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
const string DATASET = "data/test.csv";
int** matrix;
int cost;
int n;

int getCost(int* s, int** matrix, int n) {
	int cost = 0;
	int i1;
	int j1;

	for(int i = 0; i < n; i++) {
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
			printf("%4d", v[i]);
		}
		printf("\n");
	}
}

double testSwap(int *v, const int i, const int j, double cost) {
	int firstElement = v[i];
	int secondElement = v[j];
	double newcost = 0.0;

	newcost = newcost - matrix[firstElement][secondElement];
	newcost = newcost + matrix[secondElement][firstElement];

	for(int i1 = i+1; i1 < j; i1++) {
		newcost = newcost - matrix[firstElement][v[i1]];
		newcost = newcost + matrix[secondElement][v[i1]];
		newcost = newcost - matrix[v[i1]][secondElement];
		newcost = newcost + matrix[v[i1]][firstElement];
	}

	return newcost;
}

void swap(int *v, const int i, const int j) {
	int t;
	t = v[i];
	v[i] = v[j];
	v[j] = t;
}

void rotateLeft(int *v, const int start, const int n) {
	int tmp = v[start];
	for (int i = start; i < n-1; i++) {
		v[i] = v[i+1];
	}
	v[n-1] = tmp;
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

void localSearch(int* s, const int n) {
	int* sh;
	shuffle(sh, n);
	double newcost;
	bool improving = true;

	/*for(int i1 = 0; i1 < n; i1++) {
		i = sh[i1];
		rotateLeft(*s, i, n);
	}*/
	
	while(improving) {
		for(int i = 0; i < n; i++) {
			improving = false;
			for(int j = i+1; j < n; j++) {
				newcost = testSwap(s, i, j, cost);
				if(newcost < cost) {
					swap(s, i , j);
					cost = newcost;
					improving = true;
				}
			}
		}
	}
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

/*void loadData() {
	ifstream file( (INPUT_PATH + instance.file).c_str() );

	int firstColumn = 0;
	int lastColumn = d + 1;

	if(instance.hasId == true) {
		lastColumn++;
	}

    for(int row = 0; row < n; row++) {
        string line;
        getline(file, line);
        if ( !file.good() ) {
        	cout << "Error on file line reading" << endl;
            break;
        }

        stringstream iss(line);
        int j = 0;

        for(int col = firstColumn; col < lastColumn; col++) {
            string val;
            getline(iss, val, instance.delimiter);
            //data[row][j] = atof(val.c_str());
            //j++;
        }
    }
}*/

int main() {
	srand(1607);
	loadData(DATASET);

	int s0[n];
	shuffle(s0, n);

	double costS;
	double bestCost = MAX_FLOAT;
	int* bestSolution = new int[n];
	int p = 1;
	int numPert = 1000;
	int numExchanges = 0.4*n;
	double costS_;

	cost = getCost(s0, matrix, n);
	printf("f(x) = %d\n", cost);

	localSearch(s0, n);
	printf("f(x) = %d\n", cost);

	/*while(p <= numPert) {
		s = perturbation(s, numExchanges, n);
		localSearch(s);
		costS = getCost(s);
		
		if(costS < bestCost) {
			bestCost = costS;
			bestSolution = s;
		}
		p++;
	}

	//printf("g(%d) = %.15g\n", i, bestCost);*/

	return 0;
}