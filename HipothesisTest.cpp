#include <iostream>
#include <cmath>
#include <iomanip>

#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(const float nBeta, const float nMaximus, const int SIZE, const float THRESHOLD):NormalizedResidual(SIZE, THRESHOLD){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberSelectedMeasurements(0);
    setNumberNonSelectedMeasurements(0);
}

HipothesisTest::~HipothesisTest(){

    delete [] suspect_selected_measurements;
    delete [] suspect_residual_covariance_matrix;
    delete [] suspect_error_measurements;
    delete [] true_measurements;
    delete [] sensitivity_matrix_SS;
    delete [] sensitivity_matrix_ST;
    delete [] suspect_residual_measurements;
    delete [] inverse_sensitivity_matrix_ss;
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

void HipothesisTest::setNumberNonSelectedMeasurements(int numberNonSelectedMeasurements){
    number_non_selected_measurements = numberNonSelectedMeasurements;
}


//Funcoes GET ===========================================================

int HipothesisTest::getNumberSelectedMeasurements() const{
    return number_selected_measurements;
}

int HipothesisTest::getNumberNonSelectedMeasurements() const{
    return number_non_selected_measurements;
}

float HipothesisTest::setNBeta() const{
    return N_beta;
}

float HipothesisTest::setNMaximus() const{
    return N_maximus;
}

float *HipothesisTest::getNonSelectedMeasurements()const{
    return non_selected_measurements;
}

float *HipothesisTest::getSuspectSelectedMeasurements()const{
    return suspect_selected_measurements;
}

float *HipothesisTest::getSuspectResidualCovarianceMatrix()const{
    return suspect_residual_covariance_matrix;
}

float *HipothesisTest::getSuspectErrorMeasurements()const{
    return suspect_error_measurements;
}

float *HipothesisTest::getTrueMeasurements()const{
    return true_measurements;
}

float *HipothesisTest::getSensitivityMatrixSS()const{
    return sensitivity_matrix_SS;
}

float *HipothesisTest::getInverseSensitivityMatrixSS()const{
    return inverse_sensitivity_matrix_ss;
}

float *HipothesisTest::getSensitivityMatrixST()const{
    return sensitivity_matrix_ST;
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

    number_non_selected_measurements = num - number_selected_measurements;

    suspect_selected_measurements = new float[number_selected_measurements];
    non_selected_measurements = new float[number_non_selected_measurements];
    
    int count = 0;
    int count1 = 0;
    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){   
            suspect_selected_measurements[count] = i;
            count++;
        }

        else{
            non_selected_measurements[count1] = i;
            count1++; 
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

    true_measurements = new float[number_non_selected_measurements];
    
    for(int i = 0; i < number_non_selected_measurements; i++){
        j = non_selected_measurements[i];
        true_measurements[i] = temp[j];
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

void HipothesisTest::SelectSensitivityMatrixST(){

    int num = getNumberOfMeasurements();
    float *temp = getSensitivityMatrix();
    int k;

    sensitivity_matrix_ST = new float[number_selected_measurements*number_non_selected_measurements];
    //achar matriz de sensibilidade ST

    for(int i = 0; i < number_selected_measurements; i++){
        for(int j = 0; j < number_non_selected_measurements; j++){
            sensitivity_matrix_ST[i*number_non_selected_measurements + j] = 0;
        }
    }

    for(int i = 0; i < number_selected_measurements; i++){
        k = suspect_selected_measurements[i];
        for(int j = 0; j < number_selected_measurements; j++){
            for(int n = 0; n < number_non_selected_measurements; n++){
                
                if((j == k || n == k)){
                    sensitivity_matrix_ST[j*number_non_selected_measurements + n] = temp[j*num + n];
                }

            }
        }
    }

}


void HipothesisTest::CalculateInverseSensitivityMatrixSS(){

    inverse_sensitivity_matrix_ss = new float[number_selected_measurements*number_selected_measurements];

    inverse_sensitivity_matrix_ss = CalculateInverseMatrix(sensitivity_matrix_SS, number_selected_measurements);
    
}

void HipothesisTest::CalculateSuspectResidualMeasurements(){
    
    float *temp = new float[number_selected_measurements];
    
    suspect_residual_measurements = new float[number_selected_measurements];

    suspect_residual_measurements = MultiplyArray(sensitivity_matrix_SS,suspect_error_measurements,number_selected_measurements,number_selected_measurements,number_selected_measurements,1);
    temp = MultiplyArray(sensitivity_matrix_ST,true_measurements,number_selected_measurements,number_non_selected_measurements,number_non_selected_measurements,1);
    
    for(int i = 0; i < number_selected_measurements; i++){
        suspect_residual_measurements[i] = suspect_residual_measurements[i] + temp[i];
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
