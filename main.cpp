#include<iostream>
#include<cmath>
#include<iomanip>

#include"NormalizedResidual.h"

int main(){
    // const int size = 5;
    // const float thershold = 3.0;

    // const float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    // const float estimatedMeasurement[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    // const double covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
    //                                    {0,0.00003285,0,0,0},
    //                                    {0,0, 0.000006805,0,0},
    //                                    {0,0,0, 0.000006805,0},
    //                                    {0,0,0,0, 0.0000405}};

    const int size = 4;
    const float threshold = 3.0;

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

    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    NormalizedResidual LNRTest(size,threshold);
    // LNRTest.LargestNormalizedResidualTest(measurementArray,estimatedArray,jPtr,gPtr,cmPtr);
    LNRTest.setMeasurementArray(measurementArray);
    LNRTest.setEstimatedMeasurementArray(estimatedArray);
    LNRTest.setResidualCovarianceMatrix(cmPtr);

    float largestResidual;
    int position;

    for(int i=0; i < size; i++){
        LNRTest.calculateResidualArray();
        LNRTest.calculateNormalizedResidualArray();
        if(i==0){
            std::cout << "Conjunto de medicoes residuais: " << std::endl; 
            LNRTest.print(LNRTest.getResidualArray());
            std::cout << std::endl << "Conjunto de medicoes residuais normalizadas: " << std::endl; 
            LNRTest.print(LNRTest.getNormalizedArray());
            std::cout << std::endl << "Limite: " << std::setprecision(3) <<threshold << std::endl;
        }
        LNRTest.findLargestResidual(largestResidual, position);
        if (largestResidual > threshold){ 
            LNRTest.deleteError(threshold, largestResidual, position);
            LNRTest.print(LNRTest.getNormalizedArray());
        }
        else{
            std::cout << std::endl << "Conjunto de medicoes livre de erro!" << std::endl <<
            std::endl;
            break;
        }
    }

    return 0;
}

