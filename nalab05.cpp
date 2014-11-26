#include <stdio.h>
#include <math.h>
#include <iostream>
#define MAX_SIZE 100
using namespace std;

double absolute(double value);
int factorization(int n, double (&a)[MAX_SIZE][MAX_SIZE], double (&L)[MAX_SIZE][MAX_SIZE], double (&U)[MAX_SIZE][MAX_SIZE]);
void solve(int n, double (&L)[MAX_SIZE][MAX_SIZE], double (&U)[MAX_SIZE][MAX_SIZE], double x[], double b[]);
int EigenV(int n, double (&a)[MAX_SIZE][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);
 
int main()
{
    int n, MAXN, m, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
    double lambda, TOL;
 
    while (scanf("%d", &n) != EOF) {
        for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        scanf("%lf", &a[i][j]);
        scanf("%lf %d", &TOL, &MAXN);
        scanf("%d", &m);
        for (i=0; i<m; i++) {
            scanf("%lf", &lambda);
            for (j=0; j<n; j++)
            scanf("%lf", &v[j]);
            switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
                case -1:
                    printf("%12.8f is an eigenvalue.\n", lambda );
                    break;
                case 0:
                    printf("Maximum number of iterations exceeded.\n");
                    break;
                case 1:
                    printf("%12.8f\n", lambda );
                    for (k=0; k<n; k++)
                    printf("%12.8f ", v[k]);
                    printf("\n");
                    break;
            }
        }
        printf("\n");
    }
 
    return 0;
}

double absolute(double value){
    if (value<0) return (-value);
    else return value;
}

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
        U[n-1][n-1] -= L[n-1][k] * U[k][n-1];
    }
    return 0;
}

void solve(int n, double (&L)[MAX_SIZE][MAX_SIZE], double (&U)[MAX_SIZE][MAX_SIZE], double x[], double b[]){
    double y[n];
    y[0] = b[0] / L[0][0];
    for (int i=1;i<n;i++){
        y[i] = b[i];
        for (int j=0;j<i;j++){
            //cout<<"j="<<j<<endl;
            y[i] -= L[i][j] * y[j];
        }
    }
    x[n-1] = y[n-1] / U[n-1][n-1];
    for (int i=n-2;i>=0;i--){
        x[i] = y[i];
        for (int j=i+1;j<n;j++){
            x[i] -= U[i][j] *x[j];
        }
        x[i] /= U[i][i];
    }
}

int EigenV(int n, double (&a)[MAX_SIZE][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN){
    double q = *lambda;
    cout<<"q="<<q<<endl;//
    int k=1,p=0;
    double x[n];
    double _a[MAX_SIZE][MAX_SIZE];
    double L[MAX_SIZE][MAX_SIZE],U[MAX_SIZE][MAX_SIZE];
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            _a[i][j] = a[i][j];
            L[i][j] = U[i][j] = 0;
        }
    }
    /*
    int q,q1=0,q2=0;
    for (int i=0;i<n;i++){
        int t=0;
        for (int j=0;j<n;j++){
            t += v[j]*a[j][i];
        }
        q1 += v[i]*t;
        q2 += v[i]*v[i];
    }
    q = q1/q2;
    if (q == (*lambda)) return -1;*/
    
    for (int i=0;i<n;i++){
        _a[i][i] -= q;
    }

    for (int i=1;i<n;i++){
        if ( absolute(v[i]) > absolute(v[p]) ){
            p=i;
        }
    }
    for (int i=0;i<n;i++){
        v[i] /= v[p];
    }

    if (factorization(n,_a,L,U)==1){
        return -1;
    }

    cout<<"~~~~~~~~~~~"<<endl;
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout<<L[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"~~~"<<endl;
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout<<U[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"~~~~~~~~~~~"<<endl;

    //cout<<lambda<<endl;
    //cout<<(*lambda)<<"________"<<endl;
    
    while (k<=MAXN){
        solve(n,L,U,x,v);
        double niu = x[p];
        p=0;
        for (int i=1;i<n;i++){
            if ( absolute(x[i]) > absolute(x[p]) ){
                p=i;
            }
        }
        double tmp;
        double ERR = absolute(v[0] - x[0] / x[p]);
        for (int i=1;i<n;i++){
            tmp = v[i] - x[i] / x[p];
            if ( absolute(tmp) > ERR ){
                ERR = tmp;
            }
        }
        for (int i=0;i<n;i++){
            v[i] = x[i] / x[p];
        }
        if (ERR<TOL){
            niu = 1.0/niu + q;
            *lambda = niu;
            return  1; //process succeedeed
        }
        k++;
    }
    //cout<<"Maximum number of iterations exceeded"<<endl;
    return 0;
}