#include <iostream>
#include <cmath>

#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(float nBeta, float nMaximo, int SIZE){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberOfMeasurements(SIZE);

    suspect_covariance_matrix = new float[SIZE];
    suspect_residual_measurements = new float[SIZE];
    suspect_error_measuremnets = new float[SIZE];
    true_error_measurements = new float[SIZE];
    sensitivity_matrix_SS = new float[SIZE];
    inverse_sensitivity_matrix_ss = new float[SIZE];
    sensitivity_matrix_ST = new float[SIZE];
    estimated_error_measurements = new float[SIZE];
}

HipothesisTest::~HipothesisTest(){

}

//Funcoes SET ============================================================
void HipothesisTest::setNBeta(float nBeta){
    N_beta = nBeta;
}

void HipothesisTest::setNMaximus(float nMaximus){
    N_maximus = nMaximus
}

void HipothesisTest::setInverseSensitivityMatrixSS(const float *inverseSensitivityMatrixSS){
    inverse_sensitivity_matrix_ss = inverseSensitivityMatrixSS;    
}

void HipothesisTest::setSuspectCovarianceMatrix(const float *suspectCovarianceMatrix){
    suspect_covariance_matrix = suspectCovarianceMatrix;
}

void HipothesisTest::setSuspectResidualMeasurements(const float *suspectResidualMeasurements){
    suspect_residual_measurements = suspectResidualMeasurements;
}

void HipothesisTest::setSuspectErrorMeasurements(const float *suspectErrorMeasurements){
    suspect_error_measuremnets = suspectErrorMeasurements;
}

void HipothesisTest::setTrueErrorMeasurements(const float *trueErrorMeasurements){
    true_error_measurements = trueErrorMeasurements;
}

void HipothesisTest::setSensitivityMatrixSS(const float *sensitivityMatrixSS){
    sensitivity_matrix_SS = sensitivityMatrixSS;
}

void HipothesisTest::setSensitivityMatrixST(const float *sensitivityMatrixST){
    sensitivity_matrix_ST = sensitivityMatrixST;
}
void HipothesisTest::setEstimatedErrorMeasurements(const float *estimatedErrorMeasurements){
    estimated_error_measurements = estimatedErrorMeasurements;
}   

//=======================================================================

void HipothesisTest::CalculateSuspectCovarianceMatrix(const float *covarianceMatrix){

}

void HipothesisTest::CalculateSuspectErrorMeasurements(){

}

void HipothesisTest::CalculateTrueErrorMeasurements(){

}

void HipothesisTest::CalculateSuspectResidualMeasurements(){
    int num = getNumberOfMeasurements();
    for(int i = 0; i < num; i++){
        suspect_residual_measurements[i] = sensitivity_matrix_SS[i]*suspect_error_measuremnets[i];
        suspect_residual_measurements[i] += sensitivity_matrix_ST[i]*true_error_measurements[i];
    }
}

void HipothesisTest::CalculateSensitivityMatrixSS(){

}

void HipothesisTest::CalculateInverseSensitivityMatrix(const float *){

}

void HipothesisTest::CalculateSensitivityMatrixST(){

}

void HipothesisTest::CalculateEstimatedErrorMeasurements(){
    int num = getNumberOfMeasurements();
    for(int i = 0; i < num; i++){
        estimated_error_measurements[i] = inverse_sensitivity_matrix_ss[i]*suspect_residual_measurements[i];
    }
}

void HipothesisTest::CalculateNMeasurements(){
    int num = getNumberOfMeasurements();
    for(int i = 0; i < num; i++){
        N_measurements[i] = fabs(estimated_error_measurements[i]); 
        N_measurements[i] += suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]-1)*N_beta;
        N_measurements[i] /= suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
    }
}

void HipothesisTest::CalculateThersholdMeasurements(){
    int num = getNumberOfMeasurements();
    for(int i = 0; i < num; i++){
        threshold_measurements[i] = suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
        threshold_measurements[i] *= N_measurements[i];
    }
}

void HipothesisTest::SelectMeasurements(){
    int num = getNumberOfMeasurements();
    int count = 0;
    for(int i = 0; i < num; i++){
        if(estimated_error_measurements[i] > threshold_measurements[i]){
            count++;
        }
    } 

    selected_measuremnts = new float[count];

    int j = 0; 
    for(int i = 0; i < num; i++){
        if(estimated_error_measurements[i] > threshold_measurements[i]){
            selected_measurements[j] = ;
            j++
        }
    }

}
