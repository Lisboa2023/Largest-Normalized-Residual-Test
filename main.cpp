#include<iostream>
#include<cmath>
#include<iomanip>

#include"NormalizedResidual.h"
#include"HipothesisTest.h"

int main(){
    //Largest Normalized Residual Test ========================================= 
    // const int size = 5;
    // const float threshold = 3.0;

    // float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    // float estimatedMeasurement[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    // float covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
    //                                    {0,0.00003285,0,0,0},
    //                                    {0,0, 0.000006805,0,0},
    //                                    {0,0,0, 0.000006805,0},
    //                                    {0,0,0,0, 0.0000405}};

    // float *mPtr = measurementArray;
    // float *emPtr = estimatedMeasurement;
    // float *cmPtr = covarianceMatrix[0];

    // NormalizedResidual LNRTest(size,threshold);
    // LNRTest.LargestNormalizedResidualTest(mPtr,emPtr,cmPtr);
    // =============================================================================




    // Teste 2 ========================================================================
    // const int size = 4;
    // const float threshold = 3.0;

    // float measurementArray[] = {3.90, -4.05, -0.48, 2.04};
    // float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};
    // float covarianceMatrix[][size] = {{0.001, 0,0,0},
    //                                   {0,0.004,0,0},
    //                                   {0,0, 0.001,0},
    //                                   {0,0,0, 0.002}};

    // float *mPtr = measurementArray;
    // float *emPtr = estimatedArray;
    // float *cmPtr = covarianceMatrix[0];

    // NormalizedResidual LNRTest(size,threshold);
    // LNRTest.LargestNormalizedResidualTest(mPtr,emPtr,cmPtr);
    // =================================================================================





    // Hypothesis Test =================================================================
    const int size = 3;
    const float threshold = 3.0;
    int x = 0, y = 0;

    float jacobianMatrix[][2] = {{-33.33, 0},
                                {0,-20},
                                {45.8,-12.5}};
    float gainMatrix[][2] = {{38370000, -5729000},
                             {-5729000, 7812000}};
    float covarianceMatrix[][3] = {{0.000064, 0,0},
                                      {0,0.000064,0},
                                      {0,0, 0.0001}};

    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    HipothesisTest HTITest(x,y,size,threshold);

    HTITest.calculateHatMatrix(jPtr,gPtr,cmPtr,2);
    HTITest.calculateSensitivityMatrix();
    HTITest.calculateResidualCovarianceMatrix(cmPtr);
    HTITest.print(HTITest.getHatMatrix(),size,size);
    HTITest.print(HTITest.getSensitivityMatrix(),size,size);
    HTITest.print(HTITest.getResidualCovarianceMatrix(),size,size);
    // ================================================================================





    // Teste 2 =========================================================================

    // const int number_of_measurements = 4;
    // const int number_of_bus = 2;
    // const float threshold = 3.0;
    // int x = 0, y = 0;

    
    // float measurementArray[] = {3.90, -4.05, -0.48, 2.04};

    // float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};

    // float jacobianMatrix[][number_of_bus] = {{-50, -100},
    //                                          {150,-100},
    //                                          {-100, 200},
    //                                          {0, 100}};

    // float gainMatrix[][number_of_bus] = {{18125000, -18750000},
    //                                      {-18750000, 57500000}};

    // float covarianceMatrix[][number_of_measurements] = {{0.001,0,0,0},
    //                                                     {0,0.004,0,0},
    //                                                     {0,0, 0.001,0},
    //                                                     {0,0,0, 0.002}};

    // float *mPtr = measurementArray;
    // float *emPtr = estimatedArray;
    // float *jPtr = jacobianMatrix[0];
    // float *gPtr = gainMatrix[0];
    // float *cmPtr = covarianceMatrix[0];

    // HipothesisTest HTITest(x,y,number_of_measurements,threshold);


    // HTITest.calculateHatMatrix(jPtr,gPtr,cmPtr,number_of_bus);
    // std::cout << "Hat Matrix" << std::endl;
    // HTITest.print(HTITest.getHatMatrix(),number_of_measurements,number_of_measurements);

    // HTITest.calculateSensitivityMatrix();
    // std::cout << "Sensitivity Matrix" << std::endl;
    // HTITest.print(HTITest.getSensitivityMatrix(),number_of_measurements,number_of_measurements);

    // HTITest.setResidualCovarianceMatrix(cmPtr);
    // std::cout << "Residual Covariance Matrix" << std::endl;
    // HTITest.print(HTITest.getResidualCovarianceMatrix(),number_of_measurements,number_of_measurements);

    // HTITest.setMeasurementArray(mPtr);
    // HTITest.setEstimatedMeasurementArray(emPtr);
    // HTITest.calculateResidualArray();
    // std::cout << "Residual Measurements" << std::endl;
    // HTITest.print(HTITest.getResidualArray(),1,number_of_measurements);

    // HTITest.calculateNormalizedResidualArray();
    // std::cout << "Normalized Residual " << std::endl;
    // HTITest.print(HTITest.getNormalizedArray(),1,number_of_measurements);

    // HTITest.SelectSuspectMeasurements();
    // std::cout << "Suspect Selected Measurements " << std::endl;
    // HTITest.print(HTITest.getSuspectSelectedMeasurements(),1,HTITest.getNumberSelectedMeasurements());

    // std::cout << "Non Selected Measurements " << std::endl;
    // HTITest.print(HTITest.getNonSelectedMeasurements(),1,HTITest.getNumberNonSelectedMeasurements());

    // HTITest.SelectSuspectResidualCovarianceMatrix();
    // std::cout << "Suspect Residual Covariance Matrix " << std::endl;
    // HTITest.print(HTITest.getSuspectResidualCovarianceMatrix(),HTITest.getNumberSelectedMeasurements(),HTITest.getNumberSelectedMeasurements());

    // HTITest.SelectSuspectErrorMeasurements();
    // std::cout << "Suspect Error Measurements " << std::endl;
    // HTITest.print(HTITest.getSuspectErrorMeasurements(),1,HTITest.getNumberSelectedMeasurements());

    // HTITest.SelectTrueMeasurements();
    // std::cout << "True Measurements " << std::endl;
    // HTITest.print(HTITest.getTrueMeasurements(),1,HTITest.getNumberNonSelectedMeasurements());

    // HTITest.SelectSensitivityMatrixSS();
    // std::cout << "Sensitivity Matrix SS " << std::endl;
    // HTITest.print(HTITest.getSensitivityMatrixSS(),HTITest.getNumberSelectedMeasurements(),HTITest.getNumberSelectedMeasurements());

    // HTITest.SelectSensitivityMatrixST();
    // std::cout << "Sensitivity Matrix ST " << std::endl;
    // HTITest.print(HTITest.getSensitivityMatrixST(),HTITest.getNumberOfMeasurements(),HTITest.getNumberOfMeasurements());

    return 0;
}

