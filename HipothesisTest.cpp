#include <iostream>
#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(float nBeta, float nMaximo, int SIZE){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberOfMeasurements(SIZE);

    inverse_sensitivity_matrix_ss = new float[SIZE];
    suspect_residual_measurements = new float[SIZE];
    suspect_error_measuremnets = new float[SIZE];
    true_error_measurements = new float[SIZE];
    sensitivity_matrix_SS = new float[SIZE];
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

void HipothesisTest::setInverseSensitivityMatrix(const float *inverseSensitivityMatrixSS){
    inverse_sensitivity_matrix_ss = inverseSensitivityMatrixSS;    
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

void HipothesisTest::CalculateInverseSensitivityMatrix(const float *){

}

void HipothesisTest::CalculateSuspectResidualMeasurements(){

}

void HipothesisTest::CalculateSuspectErrorMeasurements(){

}

void HipothesisTest::CalculateTrueErrorMeasurements(){

}

void HipothesisTest::CalculateSensitivityMatrixSS(){

}

void HipothesisTest::CalculateSensitivityMatrixST(){

}

void HipothesisTest::CalculateEstimatedErrorMeasurements(){

}

void HipothesisTest::CalculateNMeasurements(){

}

void HipothesisTest::CalculateThersholdMeasurements(){

}
