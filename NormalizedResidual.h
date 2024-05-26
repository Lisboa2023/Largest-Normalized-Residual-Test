#ifndef NORMALIZEDRESIDUAL_H
#define NORMALIZEDRESIDUAL_H

class NormalizedResidual{
    public:
        NormalizedResidual(const int, const float);
        ~NormalizedResidual();

        //Funcoes SET ============================================================
        void setNumberOfMeasurements(const int);
        void setThershold(const float);
        void setMeasurementArray(float *);
        void setEstimatedMeasurementArray(float *);
        void setResidualArray(float *);
        void setNormalizedArray(float *);
        void setSensitivityMatrix(float *);
        void setResidualCovarianceMatrix(float *);
        //=======================================================================
        
        void calculateSensitivityMatrix(const double *);
        void calculateResidualCovarianceMatrix(const double *);
        void calculateResidualArray();
        void calculateNormalizedResidualArray();
        void findLargestResidual(float &, int &);
        void deleteError(const int, const float, const int);
        void print(const float *);

        void LargestNormalizedResidualTest(float *, float *, const double *, const double *);
        
    private:
        int size;
        float threshold;
        float *measurement;
        float *estimatedMeasurement;
        float *residualArray;
        float *normalizedArray;
        float *sensitivityMatrix;
        float *residualCovarianceMatrix;
};

#endif