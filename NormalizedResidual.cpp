#include<iostream>
#include<cmath>
#include<iomanip>
#include"NormalizedResidual.h"

using namespace std;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

NormalizedResidual::NormalizedResidual(const int SIZE, const float THRESHOLD){
    setNumberOfMeasurements(SIZE);
    setThershold(THRESHOLD);

    residualArray = new float[SIZE];
    normalizedArray = new float[SIZE];
}

NormalizedResidual::~NormalizedResidual(){
    delete [] residualArray;
    delete [] normalizedArray;
}

void NormalizedResidual::setNumberOfMeasurements(const int SIZE){
    size = SIZE;
}

void NormalizedResidual::setThershold(const float THRESHOLD){
    threshold = THRESHOLD;
}

void NormalizedResidual::calculateResidualArray(const float *measurement, const float *mEstimated){
    for(int i = 0; i < size; i++){
        residualArray[i] = measurement[i] - mEstimated[i];
    }
}

void NormalizedResidual::calculateNormalizedResidualArray(const float *r, const double *cm){
    for(int i = 0; i < size; i++){
        normalizedArray[i] = abs(r[i])/sqrt(cm[i*size + i]);
    }
}

void NormalizedResidual::findLargestResidual(float &temp, int &pos){
    temp = *normalizedArray;
    for(int i = 0; i < size; i++){
        if(normalizedArray[i] >= temp){
            temp = normalizedArray[i];
            pos = i;
        }
    }
}

void NormalizedResidual::deletError(const int threshold, const float lg, const int p){
    if(lg > threshold){
        normalizedArray[p] = 0;
    }  
}

void NormalizedResidual::LargestNormalizedResidualTest(const float *measurementArray, const float *estimatedArray, const double *covarianceMatrix){
    
    calculateResidualArray(measurementArray, estimatedArray);
    calculateNormalizedResidualArray(residualArray, covarianceMatrix);
    print(normalizedArray);

    float largestResidual;
    int position;

    for(int i=0; i<size; i++){
        findLargestResidual(largestResidual, position); 
        deletError(threshold, largestResidual, position);
    }
    print(normalizedArray);
}

void NormalizedResidual::print(const float *ptr){
    for(int i = 0; i<size; i++){
        cout << setw(10) << setprecision(2) << ptr[i] ;
    }

    cout << endl;
}