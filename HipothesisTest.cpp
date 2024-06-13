#include <iostream>
#include <cmath>

#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(const float nBeta, const float nMaximus, const int SIZE){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberOfMeasurements(SIZE);
    setNumberSelectedMeasurements(0);
}

HipothesisTest::~HipothesisTest(){

    delete [] suspected_selected_measurements;
    delete [] suspect_covariance_matrix;
    delete [] suspect_residual_measurements;
    delete [] suspect_error_measuremnets;
    delete [] true_error_measurements;
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

void SelectSuspectMeasurements(){
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
    for(int i = 0; i < num; i++){
        suspect_residual_covariance_matrix[] = temp[suspect_selected_measurements[i]*num+suspect_selected_measurements[i]];
    }
}

void HipothesisTest::SelectSuspectErrorMeasurements(){

}

void HipothesisTest::SelectTrueErrorMeasurements(){

}

void HipothesisTest::CalculateSuspectResidualMeasurements(){

    for(int i = 0; i < number_selected_measurements; i++){
        suspect_residual_measurements[i] = sensitivity_matrix_SS[i]*suspect_error_measuremnets[i];
        suspect_residual_measurements[i] += sensitivity_matrix_ST[i]*true_error_measurements[i];
    }
}

void HipothesisTest::SelectSensitivityMatrixSS(){

}

void HipothesisTest::CalculateInverseSensitivityMatrix(const float *){

}

void HipothesisTest::SelectSensitivityMatrixST(){

}

void HipothesisTest::CalculateEstimatedErrorMeasurements(){
    
    for(int i = 0; i < number_selected_measurements; i++){
        estimated_error_measurements[i] = inverse_sensitivity_matrix_ss[i]*suspect_residual_measurements[i];
    }
}

void HipothesisTest::CalculateNMeasurements(){

    for(int i = 0; i < number_selected_measurements; i++){
        N_measurements[i] = fabs(estimated_error_measurements[i]); 
        N_measurements[i] += suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]-1)*N_beta;
        N_measurements[i] /= suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
    }
}

void HipothesisTest::CalculateThersholdMeasurements(){

    for(int i = 0; i < number_selected_measurements; i++){
        threshold_measurements[i] = suspect_covariance_matrix[i]*sqrt(inverse_sensitivity_matrix_ss[i]);
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

    selected_measuremnts = new float[count];

    int j = 0; 
    for(int i = 0; i < number_selected_measurements; i++){
        if(estimated_error_measurements[i] > threshold_measurements[i]){
            selected_measurements[j] = ;
            j++
        }
    }

}
