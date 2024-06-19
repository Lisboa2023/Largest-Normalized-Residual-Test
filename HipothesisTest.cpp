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

//=======================================================================

void HipothesisTest::SelectSuspectMeasurements(){
    
    float num = getNumberOfMeasurements();
    float *temp = getNormalizedArray();
    float normalizedThreshold = getThreshold();

    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){
            number_selected_measurements++;
        }
        else{
            number_non_selected_measurements++;
        }
    }

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
    int j;

    suspect_residual_covariance_matrix = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        j = suspect_selected_measurements[i];
        suspect_residual_covariance_matrix[i] = temp[j*num + j];
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
    int j;

    sensitivity_matrix_SS = new float[number_selected_measurements];
    
    for(int i = 0; i < number_selected_measurements; i++){
        j = suspect_selected_measurements[i];
        sensitivity_matrix_SS[i] = temp[j*num + j];
    }
}

void HipothesisTest::SelectSensitivityMatrixST(){

    sensitivity_matrix_ST = new float[number_selected_measurements];
    //achar matriz de sensibilidade ST

}


void HipothesisTest::CalculateInverseSensitivityMatrixSS(){

    inverse_sensitivity_matrix_ss = new float[number_selected_measurements];
    // chama função de inversao de matriz
    calculateInverseMatrix(sensitivity_matrix_SS, getNumberOfMeasurements());
    inverse_sensitivity_matrix_ss = getInverseMatrix();

}

void HipothesisTest::CalculateSuspectResidualMeasurements(){

    suspect_residual_measurements = new float[number_selected_measurements];

    for(int i = 0; i < number_selected_measurements; i++){
        suspect_residual_measurements[i] = sensitivity_matrix_SS[i]*suspect_error_measurements[i];
        suspect_residual_measurements[i] += sensitivity_matrix_ST[i]*true_measurements[i];
    }
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

void HipothesisTest::printHTI(){

    for (int i = 0; i < number_selected_measurements; i++){
        std::cout << std::setw(10) << suspect_error_measurements[i];
    }
    std::cout << std::endl;
}