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

    measurement = new float[SIZE];
    estimatedMeasurement = new float[SIZE];
    residualArray = new float[SIZE];
    normalizedArray = new float[SIZE];
}

NormalizedResidual::~NormalizedResidual(){
    delete [] measurement;
    delete [] estimatedMeasurement;
    delete [] residualArray;
    delete [] normalizedArray;
}

void NormalizedResidual::setNumberOfMeasurements(const int SIZE){
    size = SIZE;
}

void NormalizedResidual::setThershold(const float THRESHOLD){
    threshold = THRESHOLD;
}

void NormalizedResidual::setMeasurementArray(float *measurementArray){
    measurement = measurementArray;
}

void NormalizedResidual::setEstimatedMeasurementArray(float *estimatedArray){
    estimatedMeasurement = estimatedArray;
}

void NormalizedResidual::calculateResidualArray(){
    for(int i = 0; i < size; i++){
        residualArray[i] = measurement[i] - estimatedMeasurement[i];
    }
}

void NormalizedResidual::calculateNormalizedResidualArray(const double *cm){
    for(int i = 0; i < size; i++){
        normalizedArray[i] = abs(residualArray[i])/sqrt(cm[i*size + i]);
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
        measurement[p] = 0;
        estimatedMeasurement[p] = 0;
    }  
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, const double *covarianceMatrix){

    setMeasurementArray(measurementArray);
    setEstimatedMeasurementArray(estimatedArray);

    float largestResidual;
    int position;

    for(int i=0; i<size; i++){
        calculateResidualArray();
        calculateNormalizedResidualArray(covarianceMatrix);
        findLargestResidual(largestResidual, position); 
        deletError(threshold, largestResidual, position);
        print(normalizedArray);
    }
}

void NormalizedResidual::print(const float *ptr){
    for(int i = 0; i<size; i++){
        cout << setw(10) << setprecision(3) << ptr[i] ;
    }

    cout << endl;
}
