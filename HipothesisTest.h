#ifndef HTI_H
#define HTI_H

class HipothesisTest{
    public:
        HipothesisTest();
        ~HipothesisTest();


    private:
        float N_beta;
        float N_maximo;
        float inverse_sensitivity_matrix;

};

#endif