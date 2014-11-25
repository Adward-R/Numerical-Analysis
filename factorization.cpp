#include <iostream>
#define MAX_SIZE 3
using namespace std;

int factorization(int n, double (&a)[MAX_SIZE][MAX_SIZE], double (&L)[MAX_SIZE][MAX_SIZE], double (&U)[MAX_SIZE][MAX_SIZE]){
    for (int i=0;i<n;i++){
        L[i][i] = 1;
    }
    if (a[0][0]==0){
        return 1;
        //cout<<"Factorization impossible"<<endl;
    }
    else{
        U[0][0] = a[0][0];
    }
    for (int j=1;j<n;j++){
        U[0][j] = a[0][j];
        L[j][0] = a[j][0] / U[0][0];
    }
    for (int i=1;i<n-1;i++){
        if (a[i][i]==0){
            return 1;
            //cout<<"Factorization impossible"<<endl;
        }
        else{
            U[i][i] = a[i][i];
            for (int k=0;k<i;k++){
                U[i][i] -= L[i][k] * U[k][i];
            }
        }
        for (int j=i+1;j<n;j++){
            U[i][j] = a[i][j];
            for (int k=0;k<i;k++){
                U[i][j] -= L[i][k] * U[k][j];
            }
            L[j][i] = a[j][i];
            for (int k=0;k<i;k++){
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
    U[n-1][n-1] = a[n-1][n-1];
    for (int k=0;k<n-1;k++){
        U[n-1][n-1] -= L[n][k] * U[k][n];
    }
    return 0;
}

int main(){
    double a[3][3],L[3][3],U[3][3];
    //for (int i=0;i<3;i++){
    //    cin>>a[i][0]>>a[i][1]>>a[i][2];
    //}
    int n;
    cin>>n;
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }
    
    cout<<"---"<<factorization(n,a,L,U)<<endl<<"------"<<endl;
    /*for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            cout<<L[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"-----"<<endl;
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            cout<<U[i][j]<<" ";
        }
        cout<<endl;
        }*/
   
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            double tmp=0;
            for (int k=0;k<n;k++){
                tmp += L[i][0]*U[0][j];
            }
            cout<<tmp<<" ";
        }
        cout<<endl;
    }
    return 0;
}