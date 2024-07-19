#include <iostream>
#include <cmath>
#include <iomanip>

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
    delete [] sensitivity_matrix_SS;
    delete [] suspect_residual_measurements;
    delete [] inverse_sensitivity_matrix_ss;
    delete [] estimated_error_measurements;
    delete [] N_measurements;
    delete [] threshold_measurements;
    delete [] new_suspect_selected_measurements;

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

//Funcoes GET ===========================================================

int HipothesisTest::getNumberSelectedMeasurements() const{
    return number_selected_measurements;
}

float HipothesisTest::setNBeta() const{
    return N_beta;
}

float HipothesisTest::setNMaximus() const{
    return N_maximus;
}

float *HipothesisTest::getSuspectSelectedMeasurements()const{
    return suspect_selected_measurements;
}

float *HipothesisTest::getSuspectResidualCovarianceMatrix()const{
    return suspect_residual_covariance_matrix;
}

float *HipothesisTest::getSensitivityMatrixSS()const{
    return sensitivity_matrix_SS;
}

float *HipothesisTest::getInverseSensitivityMatrixSS()const{
    return inverse_sensitivity_matrix_ss;
}

float *HipothesisTest::getSuspectResidualMeasurements()const{
    return suspect_residual_measurements;
}

//=========================================================================  


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
    }
}

void HipothesisTest::SelectSuspectResidualCovarianceMatrix(){
    int num = getNumberOfMeasurements();
    float *temp = getResidualCovarianceMatrix();
    int k;

    suspect_residual_covariance_matrix = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        for(int j = 0; j < number_selected_measurements; j++){
            if(i == j){
                k = suspect_selected_measurements[i];
                suspect_residual_covariance_matrix[i*number_selected_measurements + j] = temp[k*num + k];
            }
            else{
                suspect_residual_covariance_matrix[i*number_selected_measurements + j] = 0;
            }
        }
    }
}

void HipothesisTest::SelectSensitivityMatrixSS(){
    int num = getNumberOfMeasurements();
    float *temp = getSensitivityMatrix();
    int k;

    sensitivity_matrix_SS = new float[number_selected_measurements*number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        for(int j = 0; j < number_selected_measurements; j++){
            if(i==j){
                k = suspect_selected_measurements[i];
                sensitivity_matrix_SS[i*number_selected_measurements + j] = temp[k*num + k];
            }
            else{
                sensitivity_matrix_SS[i*number_selected_measurements + j] = 0;
            }

        }
    }
}

void HipothesisTest::CalculateInverseSensitivityMatrixSS(){

    inverse_sensitivity_matrix_ss = new float[number_selected_measurements*number_selected_measurements];

    inverse_sensitivity_matrix_ss = CalculateInverseMatrix(sensitivity_matrix_SS, number_selected_measurements);
    
}

void HipothesisTest::SelectSuspectResidualMeasurements(){
    int num = getNumberOfMeasurements();
    float *temp = getResidualArray();
    int k;
    
    suspect_residual_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        k = suspect_selected_measurements[i];
        suspect_residual_measurements[i] = temp[k];
    }
}

void HipothesisTest::CalculateEstimatedErrorMeasurements(){

    estimated_error_measurements = new float[number_selected_measurements];
    
    estimated_error_measurements = MultiplyArray(inverse_sensitivity_matrix_ss,suspect_residual_measurements,number_selected_measurements,number_selected_measurements,number_selected_measurements,1);

}

void HipothesisTest::CalculateNMeasurements(){

    N_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        N_measurements[i] = fabs(estimated_error_measurements[i]); 
        N_measurements[i] += suspect_residual_covariance_matrix[i*number_selected_measurements + i]*sqrt(inverse_sensitivity_matrix_ss[i*number_selected_measurements+i]-1)*N_beta;
        N_measurements[i] /= suspect_residual_covariance_matrix[i*number_selected_measurements+i]*sqrt(inverse_sensitivity_matrix_ss[i*number_selected_measurements+i]);
    }
}

void HipothesisTest::CalculateThersholdMeasurements(){

    threshold_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        threshold_measurements[i] = suspect_residual_covariance_matrix[i*number_selected_measurements+i]*sqrt(inverse_sensitivity_matrix_ss[i*number_selected_measurements+i]);
        threshold_measurements[i] *= N_measurements[i];
    }
}

void HipothesisTest::SelectNewSuspectMeasurements(){

    int number_suspect_measurements = 0;
    for(int i = 0; i < number_selected_measurements; i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            number_suspect_measurements++;
        }
    } 

    new_suspect_selected_measurements = new float[number_suspect_measurements];

    int j = 0; 
    for(int i = 0; i < number_selected_measurements; i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            new_suspect_selected_measurements[j] = i;
            j++;
        }
    }

    print(new_suspect_selected_measurements,1,number_suspect_measurements);

}

void HipothesisTest::HypothesisTestIdentification(float *residual_array,float *normalized_residual_array, float *sensitivity_matrix, float *residual_covariance_matrix){
    
    setSensitivityMatrix(sensitivity_matrix);
    setResidualCovarianceMatrix(residual_covariance_matrix);
    setResidualArray(residual_array);
    setNormalizedArray(normalized_residual_array);

    SelectSuspectMeasurements();
    SelectSuspectResidualCovarianceMatrix();
    SelectSuspectResidualMeasurements();
    SelectSensitivityMatrixSS();
    print(sensitivity_matrix_SS,number_selected_measurements,number_selected_measurements);

    CalculateInverseSensitivityMatrixSS();
    print(inverse_sensitivity_matrix_ss,number_selected_measurements,number_selected_measurements);
    CalculateEstimatedErrorMeasurements();
    print(estimated_error_measurements,1,number_selected_measurements);
    CalculateNMeasurements();
    print(N_measurements,1,number_selected_measurements);
    CalculateThersholdMeasurements();
    print(threshold_measurements,1,number_selected_measurements);

    SelectNewSuspectMeasurements();

}
