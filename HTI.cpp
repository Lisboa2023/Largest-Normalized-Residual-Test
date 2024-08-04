#include<iostream>
#include<cmath>
#include<iomanip>

#include"HypothesisTest.h"
#include"LargestNormalizedResidualTest.h"

int main(){

    const int number_of_measurements = 4;
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

    float *mPtr = measurementArray;
    float *emPtr = estimatedArray;
    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    HypothesisTest HTITest(n_beta,n_maximus,number_of_measurements,threshold);
    HTITest.HypothesisTestIdentification(mPtr,emPtr,jPtr,gPtr,cmPtr,number_of_bus);

    return 0;
}

