#ifndef NORMALIZEDRESIDUAL_H
#define NORMALIZEDRESIDUAL_H

class NormalizedResidual{
    public:
        NormalizedResidual(const int, const float);
        ~NormalizedResidual();

        void setNumberOfMeasurements(const int);
        void setThershold(const float);
        void calculateResidualArray(const float *, const float *);
        void calculateNormalizedResidualArray(const float *, const double *);
        void findLargestResidual(float &, int &);
        void deletError(const int, const float, const int);
        void LargestNormalizedResidualTest(const float *, const float *, const double*);
        void print(const float *);

    private:
        int size;
        float threshold;
        float *residualArray;
        float *normalizedArray;
};

#endif