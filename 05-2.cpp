#include <stdio.h>
#include <math.h>
#include <iostream>
#define MAX_SIZE 100
using namespace std;
int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);
 
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

int gauss(int n,double q,double _a[][MAX_SIZE],double x[],double b[]){
    double s[n];
    int NROW[n];
    int p;
    double a[n][n];
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (j==i)
            {
                a[i][j] = _a[i][j] - q;
            }
            else{
                a[i][j] = _a[i][j];
            }
        }
    }
    //STEP 1
    for (int i = 0; i < n; i++){
        s[i] = a[i][0];
        for (int j = 1; j < n; j++){
            if (a[i][j]>s[i]){
                s[i] = a[i][j];
            }
        }
        if (s[i]==0){
            //cout<<"No unique solution!"<<endl;
            return 1;
        }
        NROW[i] = i;
    }
    //STEP 2
    for (int i = 0; i < n-1 ; ++i)
    {
        //STEP 3
        p=i;
        for (int j = i+1; j < n; ++j)
        {
            if (absolute(a[NROW[j]][i])/s[NROW[j]] > absolute(a[NROW[p]][i])/s[NROW[p]])
            {
                p = j;
            }
        }
        //STEP 4
        if (a[NROW[p]][i] == 0)
        {
            //cout<<"No unique solution!"<<endl;
            return 1;
        }
        //STEP 5
        if (NROW[i]!=NROW[p])
        {
            int NCOPY = NROW[i];
            NROW[i] = NROW[p];
            NROW[p] = NCOPY;
        }
        //STEP 6
        for (int j = i+1; j < n; ++j)
        {
            //STEP 7
            double m = a[NROW[j]][i] / a[NROW[i]][i];
            //STEP 8
            for (int k = 0; k < n; ++k)
            {
                a[NROW[j]][k] -= m * a[NROW[i]][k];
            }
        }
    }
    //STEP 9
    if (a[NROW[n-1]][n-1]==0)
    {
        //cout<<"No unique solution!"<<endl;
        return 1;
    }
    //STEP 10
    x[n-1] = b[NROW[n-1]] / a[NROW[n-1]][n-1];
    //STEP 11
    for (int i = n-2; i >= 0; --i)
    {
        double tmp = 0;
        for (int j = i+1; j < n; ++j)
        {
            tmp += a[NROW[i]][j] * x[j];
        }
        x[i] = (b[NROW[i]] - tmp) / a[NROW[i]][i];
    }
    //STEP 12
    return 0;
}

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN){
    double q = *lambda;
    int k=1,p=0;
    double x[n];
    //double _a[MAX_SIZE][MAX_SIZE];
    //for (int i=0;i<n;i++){
    //    for (int j=0;j<n;j++){
    //        _a[i][j] = a[i][j];
    //    }
    //}
    //for (int i = 0; i < n; ++i)
    //{
    //    _a[i][i] -= q;
    //}

    for (int i=1;i<n;i++){
        if ( absolute(v[i]) > absolute(v[p]) ){
            p=i;
        }
    }
    for (int i=0;i<n;i++){
        v[i] /= v[p];
    }

/*
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
    cout<<"~~~~~~~~~~~"<<endl;*/
    
    while (k<=MAXN){
        if (gauss(n,q,a,x,v)==1){
            return -1;
        }
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