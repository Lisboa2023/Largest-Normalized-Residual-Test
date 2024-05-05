//Largest resiudal normalized

#include<iostream>
#include<cmath>
#include<iomanip>

using namespace std;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

const int tam = 5; 

void calculateResidualArray(int, float *, float *, float *);
void calculateNormalizedResidual(float *, float *, double [][tam], int);
void findLargestResidual(int, float *, float &, int &);
void deletError(int, float *, float, int);
void print(int, float *);
void print(float);
void print(int, double [][tam]);

int main(){
    const int size = 5;

    float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    float estimatedMeasurements[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    double covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
                                       {0,0.00003285,0,0,0},
                                       {0,0, 0.000006805,0,0},
                                       {0,0,0, 0.000006805,0},
                                       {0,0,0,0, 0.0000405}};

    // float residualArray[] = {0.0286, 0.0235, 0.006, -0.006, 0.007};
    float residualArray[size];
    float normalizedArray[size];

    calculateResidualArray(size, residualArray, measurementArray, estimatedMeasurements);
    calculateNormalizedResidual(normalizedArray, residualArray, covarianceMatrix, size);
    print(size, normalizedArray);

    float largestResidual;
    int position;
    const float threshold = 3.0;

    for(int i=0; i<size; i++){
        findLargestResidual(size, normalizedArray, largestResidual, position); 
        deletError(threshold, normalizedArray, largestResidual, position);
    }
    print(size,normalizedArray);

    return 0;
}

void calculateResidualArray(int size, float *r, float *measurement, float *mEstimated){
    for(int i = 0; i < size; i++){
        r[i] = measurement[i] - mEstimated[i];
    }
}

void calculateNormalizedResidual(float *rn, float *r, double cm[][tam], int size){
    for(int i = 0; i < size; i++){
        rn[i] = abs(r[i])/sqrt(cm[i][i]);
    }
}

void findLargestResidual(int size, float *array, float &temp, int &pos){
    temp = *array;
    for(int i = 0; i < size; i++){
        if(array[i] >= temp){
            temp = array[i];
            pos = i;
        }
    }
}

void deletError(int threshold, float *array, float lg, int p){
    if(lg > threshold){
        array[p] = 0;
    }  
}

void print(int size, float *ptr){
    for(int i = 0; i<size; i++){
        cout << setw(10) << setprecision(2) << *(ptr+i) ;
    }

    cout << endl;
}

void print(int size, double ptr[][tam]){
    for(int i = 0; i<size; i++){
        for(int j=0; j<size; j++){
            cout << setw(10) << setprecision(2) << ptr[i][j] << "\t";
        }

        cout << endl;
    }
}