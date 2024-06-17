#ifndef NORMALIZEDRESIDUAL_H
#define NORMALIZEDRESIDUAL_H

class NormalizedResidual{
    public:
        NormalizedResidual(const int, const float);
        ~NormalizedResidual();

        //Funcoes SET ============================================================
        void setNumberOfMeasurements(const int);
        void setThreshold(const float);

        void setHatMatrix(float *);
        void setSensitivityMatrix(float *);
        void setResidualCovarianceMatrix(float *);

        void setMeasurementArray(float *);
        void setEstimatedMeasurementArray(float *);
        void setResidualArray(float *);
        void setNormalizedArray(float *);
        //=======================================================================

        //Funcoes GET============================================================
        int getNumberOfMeasurements() const;
        float getThreshold() const;

        float *getInverseMatrix() const;
        float *getSensitivityMatrix() const;
        float *getResidualCovarianceMatrix() const;
        
        float *getResidualArray() const;
        float *getNormalizedArray() const;

        //=======================================================================
        
        void calculateInverseMatrix(float *);
        void calculateTransposedMatrix(const float *);

        void calculateHatMatrix(float *, float *, float *);
        void calculateSensitivityMatrix();
        void calculateResidualCovarianceMatrix(const float *);
        void calculateResidualArray();
        void calculateNormalizedResidualArray();
        void findLargestResidual(float &, int &);
        void deleteError(const int, const float, const int);
        void print(const float *);

        void LargestNormalizedResidualTest(float *, float *, float *, float *, float*);
        
    private:
        int size;
        float threshold;
        float *inverseMatrix;
        float *transposedMatrix;
        float *measurement;
        float *estimatedMeasurement;
        float *residualArray;
        float *normalizedArray;
        float *hatMatrix;
        float *sensitivityMatrix;
        float *residualCovarianceMatrix;
};

#endif