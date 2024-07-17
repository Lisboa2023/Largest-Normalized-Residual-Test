#include<iostream>
#include<cmath>
#include<iomanip>
#include"NormalizedResidual.h"

NormalizedResidual::NormalizedResidual(const int SIZE, const float THRESHOLD){
    setNumberOfMeasurements(SIZE);
    setThreshold(THRESHOLD);

    measurement = new float[SIZE];
    estimatedMeasurement = new float[SIZE];
    hatMatrix = new float[SIZE*SIZE];
    sensitivityMatrix = new float[SIZE*SIZE];
    residualCovarianceMatrix = new float[SIZE*SIZE];
    residualArray = new float[SIZE];
    normalizedArray = new float[SIZE];
}

NormalizedResidual::~NormalizedResidual(){

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
    number_of_measurements = SIZE;
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
    return number_of_measurements;
}

float NormalizedResidual::getThreshold() const{
    return threshold;
}

float *NormalizedResidual::getHatMatrix() const{
    return hatMatrix;
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

float *NormalizedResidual::CalculateInverseMatrix(float *matrix, const int length){
    float *inverseMatrix = new float[length*length];
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

    return inverseMatrix;

}

float *NormalizedResidual::CalculateTransposedMatrix(const float *matrix, const int ROWS, const int COLUMNS){

    float *transposedMatrix = new float[ROWS*COLUMNS];

    for(int i = 0; i < ROWS; i++){
        for(int j = 0; j < COLUMNS; j++){
            transposedMatrix[j*ROWS + i] = matrix[i*COLUMNS + j];
        }
    }

    return transposedMatrix;
}

float *NormalizedResidual::MultiplyArray(const float *array_a,const float *array_b, const int rows_a, const int columns_a, const int rows_b, const int columns_b){

    float *temp = new float[rows_a*columns_b];

    for(int i = 0; i < rows_a; i++){
        for(int j = 0; j < columns_b; j++){
            for(int k = 0; k < rows_b; k++){
                if(k == 0){
                    temp[i*columns_b + j] = array_a[i*columns_a + k]*array_b[k*columns_b + j];
                }

                else{   
                    temp[i*columns_b + j] += array_a[i*columns_a + k]*array_b[k*columns_b + j];
                }
            }
        }
    }

    return temp;
}

void NormalizedResidual::calculateHatMatrix(float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix, const int length){
    
    //hatMatrix = matriz jacobiana * matriz de ganho invertida * matriz jacobiana transposta * matriz de covariancia invertida
    float *tempInverse = CalculateInverseMatrix(gainMatrix,length);
    float *temp = MultiplyArray(jacobianMatrix,tempInverse,number_of_measurements,length,length,length);

    float *tempTransposed = CalculateTransposedMatrix(jacobianMatrix,number_of_measurements,length);
    temp = MultiplyArray(temp,tempTransposed,number_of_measurements,length,length,number_of_measurements);

    tempInverse = CalculateInverseMatrix(covarianceMatrix,number_of_measurements);
    hatMatrix = MultiplyArray(temp,tempInverse,number_of_measurements,number_of_measurements,number_of_measurements,number_of_measurements);

    //Retorna as matrizes de Ganhio e Covariancia para seus valores originais
    CalculateInverseMatrix(covarianceMatrix,number_of_measurements);
}

void NormalizedResidual::calculateSensitivityMatrix(){

    //Criando Matrix Identidade =============================================
    bool identityMatrix[number_of_measurements][number_of_measurements];
    for(int i = 0; i < number_of_measurements; i++){
        for(int j = 0; j < number_of_measurements; j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    //======================================================================

    for(int i = 0; i < number_of_measurements; i++){
        for(int j = 0; j < number_of_measurements; j++){
            sensitivityMatrix[i*number_of_measurements + j] = identityMatrix[i][j] - hatMatrix[i*number_of_measurements + j];
        }
    } 

}

void NormalizedResidual::calculateResidualCovarianceMatrix(float *covarianceMatrix){

    residualCovarianceMatrix = MultiplyArray(sensitivityMatrix,covarianceMatrix,number_of_measurements,number_of_measurements,number_of_measurements,number_of_measurements);

}

void NormalizedResidual::calculateResidualArray(){

    for(int i = 0; i < number_of_measurements; i++){
        residualArray[i] = measurement[i] - estimatedMeasurement[i];
    }

}

void NormalizedResidual::calculateNormalizedResidualArray(){

    for(int i = 0; i < number_of_measurements; i++){
        normalizedArray[i] = fabs(residualArray[i])/sqrt(residualCovarianceMatrix[i*number_of_measurements + i]);
    }

}

void NormalizedResidual::findLargestResidual(float &temp, int &pos){
    temp = *normalizedArray;
    for(int i = 0; i < number_of_measurements; i++){
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

void NormalizedResidual::print(const float *array, const int rows, const int columns){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            std::cout << std::setw(20) << array[i*columns + j];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, float *residualCovarianceMatrix){

    setMeasurementArray(measurementArray);
    setEstimatedMeasurementArray(estimatedArray);
    setResidualCovarianceMatrix(residualCovarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i < number_of_measurements; i++){
        calculateResidualArray();
        calculateNormalizedResidualArray();
        if(i==0){
            std::cout << "Conjunto de medicoes residuais: " << std::endl; 
            print(residualArray, number_of_measurements, 0);
            std::cout << std::endl << "Conjunto de medicoes residuais normalizadas: " << std::endl; 
            print(normalizedArray, number_of_measurements, 0);
            std::cout << std::endl << "Limite: " << std::setprecision(3) <<threshold << std::endl;
        }
        findLargestResidual(largestResidual, position);
        if (largestResidual > threshold){ 
            deleteError(threshold, largestResidual, position);
            print(normalizedArray, 0, number_of_measurements);
        }
        else{
            std::cout << std::endl << "Conjunto de medicoes livre de erro!" << std::endl <<
            std::endl;
            break;
        }
    }
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix,const int length){

    setMeasurementArray(measurementArray);
    setEstimatedMeasurementArray(estimatedArray);
    calculateHatMatrix(jacobianMatrix, gainMatrix, covarianceMatrix,length);
    calculateSensitivityMatrix();
    calculateResidualCovarianceMatrix(covarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i < number_of_measurements; i++){
        calculateResidualArray();
        calculateNormalizedResidualArray();
        if(i==0){
            std::cout << "Conjunto de medicoes residuais: " << std::endl; 
            print(residualArray, number_of_measurements, 0);
            std::cout << std::endl << "Conjunto de medicoes residuais normalizadas: " << std::endl; 
            print(normalizedArray, number_of_measurements, 0);
            std::cout << std::endl << "Limite: " << std::setprecision(3) <<threshold << std::endl;
        }
        findLargestResidual(largestResidual, position);
        if (largestResidual > threshold){ 
            deleteError(threshold, largestResidual, position);
            print(normalizedArray, 0, number_of_measurements);
        }
        else{
            std::cout << std::endl << "Conjunto de medicoes livre de erro!" << std::endl <<
            std::endl;
            break;
        }
    }
}
