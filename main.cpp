#include<iostream>
#include<cmath>
#include<iomanip>

#include"LargestNormalizedResidualTest.h"
#include"HypothesisTest.h"

int main(){
    //Largest Normalized Residual Test ========================================= 
    const int size = 5;
    const float threshold = 3.0;

    float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    float estimatedMeasurement[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    float covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
                                       {0,0.00003285,0,0,0},
                                       {0,0, 0.000006805,0,0},
                                       {0,0,0, 0.000006805,0},
                                       {0,0,0,0, 0.0000405}};

    float *mPtr = measurementArray;
    float *emPtr = estimatedMeasurement;
    float *cmPtr = covarianceMatrix[0];

    NormalizedResidual LNRTest(size,threshold);
    LNRTest.LargestNormalizedResidualTest(mPtr,emPtr,cmPtr);
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

    /* const int number_of_measurements = 4;
    const int number_of_bus = 2;
    const float threshold = 3.0;
    int n_beta = -2.32, n_maximus = 3.0;

    
    float measurementArray[] = {3.90, -4.05, -0.48, 2.04};

    float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};

    float jacobianMatrix[][number_of_bus] = {{-50, -100},
                                             {150,-100},
                                             {-100, 200},
                                             {0, 100}};

    float gainMatrix[][number_of_bus] = {{18125000, -18750000},
                                         {-18750000, 57500000}};

    float covarianceMatrix[][number_of_measurements] = {{0.001,0,0,0},
                                                        {0,0.004,0,0},
                                                        {0,0, 0.001,0},
                                                        {0,0,0, 0.002}};

    float sensitivityMatrix[][number_of_measurements] = {{0.257919,0.158371,0.108597,0.199095},
    {0.633484,0.669683,0.687783,-0.0723982},
    {0.108597,0.171946,0.20362,-0.126697},
    {0.39819,-0.0361991,-0.253394,0.868778}};

    float residualArray[]= {-0.092,-0.44,-0.106,-0.05};
    float normalizedArray[] = {2.9093, 6.95702, 3.35201, 1.11803};

    float *rAPtr = residualArray;
    float *nAPtr = normalizedArray;
    float *sPtr = sensitivityMatrix[0];
    float *mPtr = measurementArray;
    float *emPtr = estimatedArray;
    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    HipothesisTest HTITest(n_beta,n_maximus,number_of_measurements,threshold);
    HTITest.HypothesisTestIdentification(rAPtr,nAPtr,sPtr,cmPtr);
 */

    return 0;
}

