#include<iostream>
#include<cmath>
#include<iomanip>
#include"NormalizedResidual.h"

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

//Funcoes SET ==========================================================================================
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

void NormalizedResidual::setResidualArray(float *rArray){
    residualArray = rArray;
}

void NormalizedResidual::setNormalizedArray(float *nArray){
    normalizedArray = nArray;
}

void NormalizedResidual::setSensitivityMatrix(float *sMatrix){
    sensitivityMatrix = sMatrix;
}

void NormalizedResidual::setResidualCovarianceMatrix(float *rcMatrix){
    residualCovarianceMatrix = rcMatrix;
}
//=============================================================================================

void NormalizedResidual::calculateSensitivityMatrix(const double *hatMatrix){

    //Criando Matrix Identidade =============================================
    bool identityMatrix[size][size];
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    //======================================================================

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            sensitivityMatrix[i*size + j] = identityMatrix[i][j] - hatMatrix[i*size + j];
        }
    } 
}

void NormalizedResidual::calculateResidualCovarianceMatrix(const double *covarianceMatrix){

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                residualCovarianceMatrix[i*size + j] += sensitivityMatrix[i*size + k]*covarianceMatrix[k*size + j];
            }
        }
    }
}

void NormalizedResidual::calculateResidualArray(){
    for(int i = 0; i < size; i++){
        residualArray[i] = measurement[i] - estimatedMeasurement[i];
    }
}

void NormalizedResidual::calculateNormalizedResidualArray(){
    for(int i = 0; i < size; i++){
        normalizedArray[i] = fabs(residualArray[i])/sqrt(residualCovarianceMatrix[i*size + i]);
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

void NormalizedResidual::deleteError(const int threshold, const float lg, const int p){
    if(lg > threshold){
        residualArray[p] = 0;
        normalizedArray[p] = 0;
        measurement[p] = 0;
        estimatedMeasurement[p] = 0;
    }  
}

void NormalizedResidual::print(const float *ptr){
    for(int i = 0; i<size; i++){
        std::cout << std::setw(10) << std::setprecision(3) << ptr[i] ;
    }
    std::cout << std::endl;
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, const double *hatMatrix, const double *covarianceMatrix){

    setMeasurementArray(measurementArray);
    setEstimatedMeasurementArray(estimatedArray);
    calculateSensitivityMatrix(hatMatrix);
    calculateResidualCovarianceMatrix(covarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i<size; i++){
        calculateResidualArray();
        calculateNormalizedResidualArray();
        if(i==0){
            std::cout << "Conjunto de medicoes residuais: " << std::endl; 
            print(residualArray);
            std::cout << std::endl << "Conjunto de medicoes residuais normalizadas: " << std::endl; 
            print(normalizedArray);
            std::cout << std::endl << "Limite: " << std::setprecision(3) <<threshold << std::endl;
        }
        findLargestResidual(largestResidual, position);
        if (largestResidual > threshold){ 
            deleteError(threshold, largestResidual, position);
            print(normalizedArray);
        }
        else{
            std::cout << std::endl << "Conjunto de medicoes livre de erro!" << std::endl <<
            std::endl;
            break;
        }
    }
}
