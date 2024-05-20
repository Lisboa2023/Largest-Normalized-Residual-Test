#ifndef NORMALIZEDRESIDUAL_H
#define NORMALIZEDRESIDUAL_H

class NormalizedResidual{
    public:
        NormalizedResidual(const int, const float);
        ~NormalizedResidual();

        void setNumberOfMeasurements(const int);
        void setThershold(const float);
        void setMeasurementArray(float *);
        void setEstimatedMeasurementArray(float *);
        void calculateResidualArray();
        void calculateNormalizedResidualArray(const double *);
        void findLargestResidual(float &, int &);
        void deletError(const int, const float, const int);
        void LargestNormalizedResidualTest(float *, float *, const double*);
        void print(const float *);

    private:
        int size;
        float threshold;
        float *measurement;
        float *estimatedMeasurement;
        float *residualArray;
        float *normalizedArray;
};

#endif