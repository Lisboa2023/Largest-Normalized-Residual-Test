#include<iostream>
#include<cmath>

using namespace std;

int main(){
    const int size = 5; 
    float rn[size];
    float r[] = {4.2, 4.1, 2.3, 2.3, 1.1};
    double cm[][size] = {{0.00004637, 0,0,0,0};
                     {0,0.00003285,0,0,0};
                     {0,0, 0.000006805,0,0};
                     {0,0,0, 0.000006805,0};
                     {0,0,0,0, 0.0000405}};

    for(int i = 0; i < size; i++){
        rn[i] = abs(r[i])/sqrt(cm[i][i]);
    }

    for(int i = 0; i < size; i++){
        cout << rn[i] << endl; 
    }

    cout << endl;

    return 0;
}