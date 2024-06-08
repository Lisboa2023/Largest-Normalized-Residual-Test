#ifndef HTI_H
#define HTI_H

#include"NormalizedResidual.h"

class HipothesisTest:public NormalizedResidual{
    public:
        HipothesisTest(float, float);
        ~HipothesisTest();

        //Funcoes SET ============================================================
        void setNBeta(float);
        void setNMaximus(float);
        void setInverseSensitivityMatrix(const float *);
        void setSuspectResidualMeasurements(const float *);
        void setSuspectErrorMeasurements(const float *);
        void setTrueErrorMeasurements(const float *);
        void setSensitivityMatrixSS(const float *);
        void setSensitivityMatrixST(const float *);
        void setEstimatedErrorMeasurements(const float *);

        //=======================================================================

        void CalculateInverseSensitivityMatrix(const float *);
        void CalculateSuspectResidualMeasurements();
        void CalculateSuspectErrorMeasurements();
        void CalculateTrueErrorMeasurements();
        void CalculateSensitivityMatrixSS();
        void CalculateSensitivityMatrixST();
        void CalculateEstimatedErrorMeasurements();
        void CalculateNMeasurements();
        void CalculateThersholdMeasurements();

    private:
        float N_beta;
        float N_maximus;
        float *inverse_sensitivity_matrix_ss;
        float *suspect_residual_measurements;
        float *suspect_error_measuremnets;
        float *true_error_measurements;
        float *sensitivity_matrix_SS;
        float *sensitivity_matrix_ST;
        float *estimated_error_measurements;
        float *N_measurements;
        float *threshold_measurements;

};

#endif