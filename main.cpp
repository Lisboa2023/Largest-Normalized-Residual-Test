#include<iostream>
#include<cmath>
#include<iomanip>

#include"NormalizedResidual.h"

int main(){
    // const int size = 5;
    // const float thershold = 3.0;

    // const float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    // const float estimatedMeasurements[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    // const double covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
    //                                    {0,0.00003285,0,0,0},
    //                                    {0,0, 0.000006805,0,0},
    //                                    {0,0,0, 0.000006805,0},
    //                                    {0,0,0,0, 0.0000405}};

    const int size = 4;
    const float thershold = 3.0;

    const float measurementArray[] = {3.90, -4.05, -0.48, 2.04};
    const float estimatedMeasurements[] = {3.992, -3.61, -0.374, 2.09};
    const double covarianceMatrix[][size] = {{0.001, 0,0,0},
                                             {0,0.004,0,0},
                                             {0,0, 0.001,0},
                                             {0,0,0, 0.002}};

    const double *cmPtr = covarianceMatrix[0];

    NormalizedResidual LNRTest(size,thershold);
    LNRTest.LargestNormalizedResidualTest(measurementArray, estimatedMeasurements, cmPtr);

    return 0;
}

