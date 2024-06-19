#include<iostream>
#include<cmath>
#include<iomanip>
#include"NormalizedResidual.h"

NormalizedResidual::NormalizedResidual(const int SIZE, const float THRESHOLD){
    setNumberOfMeasurements(SIZE);
    setThreshold(THRESHOLD);

    measurement = new float[SIZE];
    estimatedMeasurement = new float[SIZE];
    sensitivityMatrix = new float[SIZE*SIZE];
    residualCovarianceMatrix = new float[SIZE*SIZE];
    residualArray = new float[SIZE];
    normalizedArray = new float[SIZE];
}

NormalizedResidual::~NormalizedResidual(){
    delete [] inverseMatrix;
    delete [] transposedMatrix;
    delete [] measurement;
    delete [] estimatedMeasurement;
    delete [] residualArray;
    delete [] normalizedArray;
    delete [] hatMatrix;
    delete [] sensitivityMatrix;
    delete [] residualCovarianceMatrix;
}

//Funcoes SET ==========================================================================================
void NormalizedResidual::setNumberOfMeasurements(const int SIZE){
    size = SIZE;
}

void NormalizedResidual::setThreshold(const float THRESHOLD){
    threshold = THRESHOLD;
}

void NormalizedResidual::setHatMatrix(float *hMatrix){
    hatMatrix = hMatrix;
}

void NormalizedResidual::setSensitivityMatrix(float *sMatrix){
    sensitivityMatrix = sMatrix;
}

void NormalizedResidual::setResidualCovarianceMatrix(float *rcMatrix){
    residualCovarianceMatrix = rcMatrix;
}

void NormalizedResidual::setMeasurementArray(float *measurementArray){
    measurement = measurementArray;
}

void NormalizedResidual::setEstimatedMeasurementArray(float *estimatedArray){
    estimatedMeasurement = estimatedArray;
}

void NormalizedResidual::setResidualArray(float *rArray){
    residualArray = rArray;
}

void NormalizedResidual::setNormalizedArray(float *nArray){
    normalizedArray = nArray;
}


//=============================================================================================


//Funcoes GET============================================================
int NormalizedResidual::getNumberOfMeasurements() const{
    return size;
}

float NormalizedResidual::getThreshold() const{
    return threshold;
}

float *NormalizedResidual::getInverseMatrix() const{
    return inverseMatrix;
}

float *NormalizedResidual::getSensitivityMatrix() const{
    return sensitivityMatrix;
}

float *NormalizedResidual::getResidualCovarianceMatrix() const{
    return residualCovarianceMatrix;
}

float *NormalizedResidual::getResidualArray() const{
    return residualArray;
}

float *NormalizedResidual::getNormalizedArray() const{
    return normalizedArray;
}

//=======================================================================

void NormalizedResidual::calculateInverseMatrix(float *matrix, const int length){

    inverseMatrix = matrix;

    //permutacao de linhas
    float temp;
    for(int i = 0; i < length; i++)
    { 
        for(int j = 0; j < length; j++)
        {
            if(i==j && inverseMatrix[i*length + j] == 0)
            {   
                for(int p = i+1; p < length; p++)
                {   
                    if(inverseMatrix[p*length + j] != 0)
                    {
                        for(int n = 0; n < length; n++)
                        {
                            temp = inverseMatrix[i*length + n];
                            inverseMatrix[i*length + n] = inverseMatrix[p*length + n];
                            inverseMatrix[p*length + n] = temp;
                        }
                        break;
                    }
                }
            }
        }
    }

    //algoritmo de eliminacao de gauss-jordan
    float identityMatrix[length][length];

    
    //Criando matriz identidade
    for(int i = 0; i < length; i++){
        for(int j = 0;j < length; j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    
    for (int i = 0; i < length; i++)
	{
		// Dividindo a linha atual pelo elemento diagonal correspondente
		float pivot = inverseMatrix[i*length + i];
		for ( int j = 0; j < length; j++)
		{
			inverseMatrix[i*length + j] /= pivot;
			identityMatrix[i][j] /= pivot; //As mesmas operações são feitas na matriz identidade para se obter a inversa
		}

		// Reduzindo as outras linhas
		for (int j = 0; j < length; j++)
		{
			if (j != i)
			{
				float a = inverseMatrix[j*length + i];
				for (int k = 0; k < length; k++)
				{
					inverseMatrix[j*length + k] -= a * inverseMatrix[i*length + k]; //Colocando 0 abaixo e acima dos pivos
					identityMatrix[j][k] -= a * identityMatrix[i][k];//Repetindo operações na matriz identidade
				}
			}
		}
	}

    for(int i = 0; i < length; i++){
        for(int j = 0; j < length; j++){
            inverseMatrix[i*length + j] = identityMatrix[i][j];
        }
    }

}

void NormalizedResidual::calculateTransposedMatrix(const float *matrix, const int ROWS, const int COLUMNS){
    
    int rows = ROWS;
    int columns = COLUMNS;

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            transposedMatrix[i*rows + j] = matrix[j*rows + i];
            std::cout << std::setw(10) << transposedMatrix[i*rows+j];
        }
    }
}

void NormalizedResidual::calculateHatMatrix(float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix){
    //hatMatrix = matriz jacobiana * matriz de ganho invertida * matriz jacobiana transposta * matriz de covariancia invertida
    hatMatrix = new float[size];
    //matriz jacobiana * matriz de ganho invertida
    calculateInverseMatrix(gainMatrix, size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                hatMatrix[i*size + j] = 0;
                hatMatrix[i*size + j] += jacobianMatrix[i*size + k]*inverseMatrix[k*size + j];
            }
        }
    }

    // * matriz jacobiana transposta
    calculateTransposedMatrix(jacobianMatrix, size,size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                hatMatrix[i*size + j] += hatMatrix[i*size + k]*transposedMatrix[k*size + j];
            }
        }
    }

    // * matriz de covariancia invertida
    calculateInverseMatrix(covarianceMatrix,size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                hatMatrix[i*size + j] += covarianceMatrix[i*size + k]*inverseMatrix[k*size + j];
            }
        }
    }
    
}

void NormalizedResidual::calculateSensitivityMatrix(){

    //Criando Matrix Identidade =============================================
    bool identityMatrix[size][size];
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    //======================================================================

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            sensitivityMatrix[i*size + j] = identityMatrix[i][j] - hatMatrix[i*size + j];
        }
    } 
}

void NormalizedResidual::calculateResidualCovarianceMatrix(const float *covarianceMatrix){

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                residualCovarianceMatrix[i*size + j] = 0;
                residualCovarianceMatrix[i*size + j] += sensitivityMatrix[i*size + k]*covarianceMatrix[k*size + j];
            }
        }
    }
}

void NormalizedResidual::calculateResidualArray(){
    for(int i = 0; i < size; i++){
        residualArray[i] = measurement[i] - estimatedMeasurement[i];
    }
}

void NormalizedResidual::calculateNormalizedResidualArray(){
    for(int i = 0; i < size; i++){
        normalizedArray[i] = fabs(residualArray[i])/sqrt(residualCovarianceMatrix[i*size + i]);
    }
}

void NormalizedResidual::findLargestResidual(float &temp, int &pos){
    temp = *normalizedArray;
    for(int i = 0; i < size; i++){
        if(normalizedArray[i] >= temp){
            temp = normalizedArray[i];
            pos = i;
        }
    }
}

void NormalizedResidual::deleteError(const int threshold, const float lg, const int p){
    if(lg > threshold){
        residualArray[p] = 0;
        normalizedArray[p] = 0;
        measurement[p] = 0;
        estimatedMeasurement[p] = 0;
    }  
}

void NormalizedResidual::print(const float *ptr){
    for(int i = 0; i < size; i++){
        std::cout << std::setw(10) << std::setprecision(3) << ptr[i] ;
    }
    std::cout << std::endl;
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix){

    setMeasurementArray(measurementArray);
    setEstimatedMeasurementArray(estimatedArray);
    setResidualCovarianceMatrix(covarianceMatrix);
    // calculateHatMatrix(jacobianMatrix, gainMatrix, covarianceMatrix);
    // calculateSensitivityMatrix();
    // calculateResidualCovarianceMatrix(covarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i < size; i++){
        calculateResidualArray();
        calculateNormalizedResidualArray();
        if(i==0){
            std::cout << "Conjunto de medicoes residuais: " << std::endl; 
            print(residualArray);
            std::cout << std::endl << "Conjunto de medicoes residuais normalizadas: " << std::endl; 
            print(normalizedArray);
            std::cout << std::endl << "Limite: " << std::setprecision(3) <<threshold << std::endl;
        }
        findLargestResidual(largestResidual, position);
        if (largestResidual > threshold){ 
            deleteError(threshold, largestResidual, position);
            print(normalizedArray);
        }
        else{
            std::cout << std::endl << "Conjunto de medicoes livre de erro!" << std::endl <<
            std::endl;
            break;
        }
    }
}
