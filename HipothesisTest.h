#ifndef HTI_H
#define HTI_H

#include"NormalizedResidual.h"

class HipothesisTest:public NormalizedResidual{
    public:
        HipothesisTest(const float, const float, const int, const float);
        ~HipothesisTest();

        //Funcoes SET ============================================================
        void setNBeta(float);
        void setNMaximus(float);
        void setNumberSelectedMeasurements(int);
        void setNumberNonSelectedMeasurements(int);

        //=======================================================================       

        void SelectSuspectMeasurements();

        void SelectSuspectResidualCovarianceMatrix();
        void SelectSuspectErrorMeasurements();
        void SelectTrueMeasurements();
        void SelectSensitivityMatrixSS();
        void SelectSensitivityMatrixST();
        
        void CalculateInverseSensitivityMatrixSS();
        void CalculateSuspectResidualMeasurements();
        void CalculateEstimatedErrorMeasurements();
        void CalculateNMeasurements();
        void CalculateThersholdMeasurements();

        void SelectMeasurements();
        void printHTI();

    private:
        int number_selected_measurements;
        int number_non_selected_measurements;
        float N_beta;
        float N_maximus;
        float *non_selected_measurements;
        float *suspect_selected_measurements;
        float *suspect_residual_covariance_matrix;
        float *suspect_residual_measurements;
        float *suspect_error_measurements;
        float *true_measurements;
        float *sensitivity_matrix_SS;
        float *inverse_sensitivity_matrix_ss;
        float *sensitivity_matrix_ST;
        float *estimated_error_measurements;
        float *N_measurements;
        float *threshold_measurements;
        float *selected_measurements;

};

#endif