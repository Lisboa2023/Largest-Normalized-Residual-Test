//Largest resiudal normalized

#include<iostream>

using namespace std;

void detectError(int, float[], float &, int &);
void delError(int, int, float[], float*, float, int p);
void print(int, float[]);

int main(){
    const int size = 5;
    const float threshold = 3.0;
    int count = 0;
    float array[] = {4.2, 4.1, 2.3, 2.3, 1.1};
    float largestResidual;
    int p;

    detectError(size, array, largestResidual, p);
    
    int tam = size - count;
    float *ptr = new float[tam];
    
    delError(size, threshold, array, ptr, largestResidual, p);
    print(tam,ptr);

    return 0;
}

void detectError(int size, float array[], float &temp, int &position){
    temp = array[0];
    for(int i = 0; i < size; i++){
        if(array[i] >= temp){
            temp = array[i];
            position = i;
        }
    }
}

void delError(int size, int threshold, float array[], float *ptr, float lg, int p){
    int j = 0;
    for (int i = 0; i < size; i++){
        if(lg > threshold){
            if(i != p)
            ptr[j] = array[i];
            j++;
        }
    }
}

void print(int tam, float ptr[]){
    for(int i = 0; i<tam; i++){
        cout << ptr[i] << endl;
    }
}