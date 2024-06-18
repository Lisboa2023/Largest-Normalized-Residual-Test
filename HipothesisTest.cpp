#include <iostream>
#include <cmath>

#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(const float nBeta, const float nMaximus, const int SIZE, const float THRESHOLD):NormalizedResidual(SIZE, THRESHOLD){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberSelectedMeasurements(0);
}

HipothesisTest::~HipothesisTest(){

    delete [] suspect_selected_measurements;
    delete [] suspect_residual_covariance_matrix;
    delete [] suspect_residual_measurements;
    delete [] suspect_error_measurements;
    delete [] true_measurements;
    delete [] sensitivity_matrix_SS;
    delete [] inverse_sensitivity_matrix_ss;
    delete [] sensitivity_matrix_ST;
    delete [] estimated_error_measurements;
    delete [] N_measurements;
    delete [] threshold_measurements;
    delete [] selected_measurements;

}

//Funcoes SET ============================================================
void HipothesisTest::setNBeta(float nBeta){
    N_beta = nBeta;
}

void HipothesisTest::setNMaximus(float nMaximus){
    N_maximus = nMaximus;
}

void HipothesisTest::setNumberSelectedMeasurements(int numberSelectedMeasurements){
    number_selected_measurements = numberSelectedMeasurements;
}

void HipothesisTest::setInverseSensitivityMatrixSS(float *inverseSensitivityMatrixSS){
    inverse_sensitivity_matrix_ss = inverseSensitivityMatrixSS;    
}

void HipothesisTest::setSuspectCovarianceMatrix(float *suspectCovarianceMatrix){
    suspect_residual_covariance_matrix = suspectCovarianceMatrix;
}

void HipothesisTest::setSuspectResidualMeasurements(float *suspectResidualMeasurements){
    suspect_residual_measurements = suspectResidualMeasurements;
}

void HipothesisTest::setSuspectErrorMeasurements(float *suspectErrorMeasurements){
    suspect_error_measurements = suspectErrorMeasurements;
}

void HipothesisTest::setTrueMeasurements(float *trueMeasurements){
    true_measurements = trueMeasurements;
}

void HipothesisTest::setSensitivityMatrixSS(float *sensitivityMatrixSS){
    sensitivity_matrix_SS = sensitivityMatrixSS;
}

void HipothesisTest::setSensitivityMatrixST(float *sensitivityMatrixST){
    sensitivity_matrix_ST = sensitivityMatrixST;
}
void HipothesisTest::setEstimatedErrorMeasurements(float *estimatedErrorMeasurements){
    estimated_error_measurements = estimatedErrorMeasurements;
}   

//=======================================================================

void HipothesisTest::SelectSuspectMeasurements(){
    
    float num = getNumberOfMeasurements();
    float *temp = getNormalizedArray();
    float normalizedThreshold = getThreshold();

    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){
            number_selected_measurements++;
        }
    }

    suspect_selected_measurements = new float[number_selected_measurements];
    
    int count = 0;
    int count1 = 0;
    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){   
            suspect_selected_measurements[count] = i;
            count++;
        }

        else{
            // true_selected_measurements[count1] = i;
            // count1++; 
        } 
    }
}

void HipothesisTest::SelectSuspectResidualCovarianceMatrix(){
    float *temp = getResidualCovarianceMatrix();
    int j;

    suspect_residual_covariance_matrix = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        j = suspect_selected_measurements[i];
        suspect_residual_covariance_matrix[i] = temp[j*number_selected_measurements + j];
    }
}

void HipothesisTest::SelectSuspectErrorMeasurements(){
    float *temp = getResidualArray();
    int j;
    suspect_error_measurements = new float[number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        j = suspect_selected_measurements[i];
        suspect_error_measurements[i] = temp[j];
    }
}

void HipothesisTest::SelectTrueMeasurements(){
    float *temp = getResidualArray();
    int j;

    true_measurements = new float[number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        j = true_selected_measurements[i];
        true_measurements[i] = temp[j];
    }
}

void HipothesisTest::CalculateSuspectResidualMeasurements(){

    suspect_residual_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        suspect_residual_measurements[i] = sensitivity_matrix_SS[i]*suspect_error_measurements[i];
        suspect_residual_measurements[i] += sensitivity_matrix_ST[i]*true_measurements[i];
    }
}


void HipothesisTest::SelectSensitivityMatrixSS(){
    float *temp = getSensitivityMatrix();
    int j;

    sensitivity_matrix_SS = new float[number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        j = suspect_selected_measurements[i];
        sensitivity_matrix_SS[i] = temp[j*number_selected_measurements + j];
    }
}

void HipothesisTest::CalculateInverseSensitivityMatrixSS(){

    inverse_sensitivity_matrix_ss = new float[number_selected_measurements];
    // chama função de inversao de matriz
    calculateInverseMatrix(sensitivity_matrix_SS);
    inverse_sensitivity_matrix_ss = getInverseMatrix();

}

void HipothesisTest::SelectSensitivityMatrixST(){

    sensitivity_matrix_ST = new float[number_selected_measurements];
    //achar matriz de sensibilidade ST

}

void HipothesisTest::CalculateEstimatedErrorMeasurements(){

    estimated_error_measurements = new float[number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        estimated_error_measurements[i] = inverse_sensitivity_matrix_ss[i]*suspect_residual_measurements[i];
    }
}

void HipothesisTest::CalculateNMeasurements(){

    N_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        N_measurements[i] = fabs(estimated_error_measurements[i]); 
        N_measurements[i] += suspect_residual_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]-1)*N_beta;
        N_measurements[i] /= suspect_residual_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
    }
}

void HipothesisTest::CalculateThersholdMeasurements(){

    threshold_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        threshold_measurements[i] = suspect_residual_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
        threshold_measurements[i] *= N_measurements[i];
    }
}

void HipothesisTest::SelectMeasurements(){

    int count = 0;
    for(int i = 0; i < number_selected_measurements; i++){
        if(estimated_error_measurements[i] > threshold_measurements[i]){
            count++;
        }
    } 

    selected_measurements = new float[count];

    int j = 0; 
    for(int i = 0; i < number_selected_measurements; i++){
        if(estimated_error_measurements[i] > threshold_measurements[i]){
            selected_measurements[j] = i;
            j++;
        }
    }

}
