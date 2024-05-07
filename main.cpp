//Largest resiudal normalized

#include<iostream>
#include<cmath>
#include<iomanip>

using namespace std;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

void calculateResidualArray(float *, const float *, const float *);
void calculateNormalizedResidual(float *, const float *, const double *);
void findLargestResidual(const float *, float &, int &);
void deletError(const int, float *, const float, const int);
void print(const float *);
const int size = 5;
int main(){
    const int size = 5;

    const float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    const float estimatedMeasurements[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    const double covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
                                       {0,0.00003285,0,0,0},
                                       {0,0, 0.000006805,0,0},
                                       {0,0,0, 0.000006805,0},
                                       {0,0,0,0, 0.0000405}};

    const double *cmPtr = covarianceMatrix[0];

    float residualArray[size];
    float normalizedArray[size];

    calculateResidualArray(residualArray, measurementArray, estimatedMeasurements);
    calculateNormalizedResidual(normalizedArray, residualArray, cmPtr);
    print(normalizedArray);

    float largestResidual;
    int position;
    const float threshold = 3.0;

    for(int i=0; i<size; i++){
        findLargestResidual(normalizedArray, largestResidual, position); 
        deletError(threshold, normalizedArray, largestResidual, position);
    }
    print(normalizedArray);

    return 0;
}

void calculateResidualArray(float *r, const float *measurement, const float *mEstimated){
    for(int i = 0; i < size; i++){
        r[i] = measurement[i] - mEstimated[i];
    }
}

void calculateNormalizedResidual(float *rn, const float *r, const double *cm){
    for(int i = 0; i < size; i++){
        rn[i] = abs(r[i])/sqrt(cm[i*size + i]);
    }
}

void findLargestResidual(const float *array, float &temp, int &pos){
    temp = *array;
    for(int i = 0; i < size; i++){
        if(array[i] >= temp){
            temp = array[i];
            pos = i;
        }
    }
}

void deletError(const int threshold, float *array, const float lg, const int p){
    if(lg > threshold){
        array[p] = 0;
    }  
}

void print(const float *ptr){
    for(int i = 0; i<size; i++){
        cout << setw(10) << setprecision(2) << ptr[i] ;
    }

    cout << endl;
}