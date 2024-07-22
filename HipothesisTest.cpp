#include <iostream>
#include <cmath>
#include <iomanip>

#include "NormalizedResidual.h"
#include "HipothesisTest.h"

HipothesisTest::HipothesisTest(const float nBeta, const float nMaximus, const int SIZE, const float THRESHOLD):NormalizedResidual(SIZE, THRESHOLD){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberSelectedMeasurements(0);
    new_number_suspect_measurements = 0;
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
   
    for(int i = 0; i < number_selected_measurements; i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            new_number_suspect_measurements++;
        }
    } 

    float *new_suspect_selected_measurements = new float[new_number_suspect_measurements];

    int j = 0; 
    for(int i = 0; i < number_selected_measurements; i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            new_suspect_selected_measurements[j] = i;
            j++;
        }
    }

    print(new_suspect_selected_measurements,1,new_number_suspect_measurements);

    float *new_suspect_residual_measurements = new float[new_number_suspect_measurements];
    float *new_suspect_residual_covariance_matrix = new float[new_number_suspect_measurements*new_number_suspect_measurements];
    float *new_sensitivity_matrix_ss = new float[new_number_suspect_measurements*new_number_suspect_measurements];

    int k;
    if(number_selected_measurements > new_number_suspect_measurements){
        
        //Novas medicoes residuais
        for(int i = 0; i < new_number_suspect_measurements; i++){
            k = new_suspect_selected_measurements[i];
            new_suspect_residual_measurements[i] = suspect_residual_measurements[k];
        }

        delete [] suspect_residual_measurements;

        suspect_residual_measurements = new float[new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            suspect_residual_measurements[i] = new_suspect_residual_measurements[i];
        }

        //Nova matrix de covariancia residual
        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                if(i == j){
                    k = new_suspect_selected_measurements[i];
                    new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = suspect_residual_covariance_matrix[k*number_selected_measurements + k];
                }

                else{
                    new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = 0;
                }
            }
        }

        delete [] suspect_residual_covariance_matrix;

        suspect_residual_covariance_matrix = new float [new_number_suspect_measurements*new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j];
            }
        }

        //Nova matrix de sensibilidade residual
        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                if(i==j){
                    k = new_suspect_selected_measurements[i];
                    new_sensitivity_matrix_ss[i*new_number_suspect_measurements + j] = sensitivity_matrix_SS[k*number_selected_measurements + k];
                }

                else{
                    new_sensitivity_matrix_ss[i*new_number_suspect_measurements + j] = 0;
                }
            }
        }

        delete [] sensitivity_matrix_SS;

        sensitivity_matrix_SS = new float[new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                sensitivity_matrix_SS[i*new_number_suspect_measurements + j] = new_sensitivity_matrix_ss[i*number_selected_measurements + j];
            }
        }

        //Novo numero de mediçoes
        number_selected_measurements = new_number_suspect_measurements;

        delete [] new_suspect_selected_measurements;
        delete [] new_suspect_residual_measurements;
        delete [] new_suspect_residual_covariance_matrix;
        delete [] new_sensitivity_matrix_ss;

    }

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
    std::cout << "Sensitivity Matrix SS" << std::endl;
    print(sensitivity_matrix_SS,number_selected_measurements,number_selected_measurements);

    while (number_selected_measurements > new_number_suspect_measurements)
    {
    
        CalculateInverseSensitivityMatrixSS();
        std::cout << "Inverse Sensitivity Matrix SS" << std::endl;
        print(inverse_sensitivity_matrix_ss,number_selected_measurements,number_selected_measurements);
        CalculateEstimatedErrorMeasurements();
        std::cout << "Estimated Error" << std::endl;
        print(estimated_error_measurements,1,number_selected_measurements);
        CalculateNMeasurements();
        std::cout << "N measurements" << std::endl;
        print(N_measurements,1,number_selected_measurements);
        CalculateThersholdMeasurements();
        std::cout << "Threshold" << std::endl;    
        print(threshold_measurements,1,number_selected_measurements);
        std::cout << "New Suspect Measurements" << std::endl;
        SelectNewSuspectMeasurements();
        
    }
    
}
