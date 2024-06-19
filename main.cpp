#include<iostream>
#include<cmath>
#include<iomanip>

#include"NormalizedResidual.h"
#include"HipothesisTest.h"

int main(){
    // const int size = 5;
    // const float thershold = 3.0;

    // float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    // float estimatedMeasurement[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    // float covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
    //                                    {0,0.00003285,0,0,0},
    //                                    {0,0, 0.000006805,0,0},
    //                                    {0,0,0, 0.000006805,0},
    //                                    {0,0,0,0, 0.0000405}};

    const int size = 4;
    const float threshold = 3.0;
    int x =0, y = 0;

    float measurementArray[] = {3.90, -4.05, -0.48, 2.04};
    float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};
    float jacobianMatrix[][2] = {{-50, -100},
                                {150,-100},
                                {-100,200},
                                {0, -100}};
    float gainMatrix[][2] = {{1812500, -18750000},
                             {-18750000, 57500000}};
    float covarianceMatrix[][size] = {{0.001, 0,0,0},
                                      {0,0.004,0,0},
                                      {0,0, 0.001,0},
                                      {0,0,0, 0.002}};

    float sensitivityMatrix[][3] = {{0.4916,-0.2236,0.3577},
                                    {-0.2236,0.1017,-0.1627},
                                    {0.5589,-0.2542,0.4068}}; 

    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    HipothesisTest HTITest(x,y,size,threshold);
    HTITest.setSensitivityMatrix(jPtr);
    HTITest.calculateTransposedMatrix(HTITest.getSensitivityMatrix(),4,2);


    return 0;
}

