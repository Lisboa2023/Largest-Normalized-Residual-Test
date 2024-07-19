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

        //=======================================================================     

        //Funcoes GET ===========================================================
        int getNumberSelectedMeasurements() const;
        int getNumberNonSelectedMeasurements() const;
        float setNBeta() const;
        float setNMaximus() const;

        float *getSuspectSelectedMeasurements()const;
        float *getSuspectResidualCovarianceMatrix()const;
        float *getSensitivityMatrixSS()const;
        float *getInverseSensitivityMatrixSS()const;
        float *getSuspectResidualMeasurements()const;

        //=========================================================================  

        void SelectSuspectMeasurements();

        void SelectSuspectResidualCovarianceMatrix();
        void SelectSensitivityMatrixSS();
        void SelectSuspectResidualMeasurements();
        
        void CalculateInverseSensitivityMatrixSS();
        void CalculateEstimatedErrorMeasurements();
        void CalculateNMeasurements();
        void CalculateThersholdMeasurements();

        void SelectNewSuspectMeasurements();

        void HypothesisTestIdentification(float*,float*,float*,float*);

    private:
        int number_selected_measurements;
        float N_beta;
        float N_maximus;
        float *suspect_selected_measurements;
        float *suspect_residual_covariance_matrix;
        float *sensitivity_matrix_SS;
        float *inverse_sensitivity_matrix_ss;
        float *suspect_residual_measurements;
        float *estimated_error_measurements;
        float *N_measurements;
        float *threshold_measurements;
        float *new_suspect_selected_measurements;

};

#endif