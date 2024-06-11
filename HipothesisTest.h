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

        void setSuspectCovarianceMatrix(const float *);

        void setSuspectErrorMeasurements(const float *);
        void setTrueErrorMeasurements(const float *);
        void setSuspectResidualMeasurements(const float *);

        void setSensitivityMatrixSS(const float *);
        void setInverseSensitivityMatrixSS(const float *);
        void setSensitivityMatrixST(const float *);

        void setEstimatedErrorMeasurements(const float *);

        //=======================================================================       

        void CalculateSuspectCovarianceMatrix(const float *);

        void CalculateSuspectErrorMeasurements();
        void CalculateTrueErrorMeasurements();
        void CalculateSuspectResidualMeasurements();
        
        void CalculateSensitivityMatrixSS();
        void CalculateInverseSensitivityMatrixSS(const float *);
        void CalculateSensitivityMatrixST();
        
        void CalculateEstimatedErrorMeasurements();
        void CalculateNMeasurements();
        void CalculateThersholdMeasurements();

        void SelectMeasurements();

    private:
        float N_beta;
        float N_maximus;
        float *suspect_covariance_matrix;
        float *suspect_residual_measurements;
        float *suspect_error_measuremnets;
        float *true_error_measurements;
        float *sensitivity_matrix_SS;
        float *inverse_sensitivity_matrix_ss;
        float *sensitivity_matrix_ST;
        float *estimated_error_measurements;
        float *N_measurements;
        float *threshold_measurements;
        float *selected_measurements;

};

#endif