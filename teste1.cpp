#include<iostream>
#include<cmath>

using namespace std;

int main(){
    int size = 5;
    float mrVector[] = {4.2, 4.1, 2.3, 2.3, 1.1};
    float covarianceMatrix[size][size];

    for(int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            if (i==j){
                covarianceMatrix[i][j] = 4; 
            }
            else{
                covarianceMatrix[i][j] = 0;
            }
        }
    }

    for (int j = 0; j < size; j++){
        cout << mrVector[j] << "\t"; 
    }

    cout << endl << endl;

    for(int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
                cout << covarianceMatrix[i][j] << "\t"; 
        }
        cout << endl;
    }

    cout << endl << endl;

    float normalizedVector[size];

    for(int i = 0; i < size; i++){
        normalizedVector[i] = mrVector[i]/sqrt(covarianceMatrix[i][i]); 
    }

    for (int j = 0; j < size; j++){
        cout << normalizedVector[j] << "\t"; 
    }

    cout << endl << endl;
    return 0;
}