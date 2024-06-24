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

        float *getHatMatrix() const;
        float *getSensitivityMatrix() const;
        float *getResidualCovarianceMatrix() const;
        
        float *getResidualArray() const;
        float *getNormalizedArray() const;

        //=======================================================================
        
        float *CalculateInverseMatrix(float *, const int);
        float *CalculateTransposedMatrix(const float *, const int, const int);
        float *MultiplyArray(const float *, const float *, const int, const int, const int, const int);

        void calculateHatMatrix(float *, float *,float *, const int);
        void calculateSensitivityMatrix();
        void calculateResidualCovarianceMatrix(float *);

        void calculateResidualArray();
        void calculateNormalizedResidualArray();

        void findLargestResidual(float &, int &);
        void deleteError(const int, const float, const int);
        void print(const float*, const int, const int);

        void LargestNormalizedResidualTest(float *, float *, float *);
        void LargestNormalizedResidualTest(float *, float *, float *, float *, float*, const int);
        
    private:
        int number_of_measurements;
        float threshold;

        float *hatMatrix;
        float *sensitivityMatrix;
        float *residualCovarianceMatrix;

        float *measurement;
        float *estimatedMeasurement;
        float *residualArray;
        float *normalizedArray;

};

#endif