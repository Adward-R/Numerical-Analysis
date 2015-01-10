#include <stdio.h>
#include <math.h>
#include <iostream>
#define MAX_m 200
#define MAX_n 5
using namespace std;
 
double f1(double x)
{
    return sin(x);
}
 
double f2(double x)
{
    return exp(x);
}
 
int OPA( double (*f)(double t), int m, double x[], double w[],
double c[], double *eps );
 
void print_results( int n, double c[], double eps)
{  
    int i;
 
    printf("%d\n", n);
    for (i=0; i<=n; i++)
       printf("%8.4e ", c[i]);
    printf("\n");
    printf("error = %6.2e\n", eps);
    printf("\n");
}
 
int main()
{
    int m, i, n;
    double x[MAX_m], w[MAX_m], c[MAX_n+1], eps;
 
    m = 90;
    for (i=0; i<m; i++) {
       x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
       w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f1, m, x, w, c, &eps);
    print_results(n, c, eps);
 
    m = 200;
    for (i=0; i<m; i++) {
       x[i] = 0.01*(double)i;
       w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f2, m, x, w, c, &eps);
    print_results(n, c, eps);
 
    return 0;
}


int OPA( double (*f)(double t), int m, double x[], double w[],double c[], double *eps ){
    double err;
    double fai0[MAX_n+1],fai1[MAX_n+1],fai2[MAX_n+1];
    double y[m];
    //double B[MAX_n+1],C[MAX_n+1];
    double B,C;
    double a[MAX_n+1];
    int n=0;

    for (int i=0;i<m;i++) {
        y[i] = f(x[i]);
    }

    //STEP 1
    //fai0(x)=1
    for (int i=0;i<=MAX_n;i++){
        fai0[i] = 0;
        fai1[i] = 0;
        fai2[i] = 0;
    }
    fai0[0] = 1;
    //a[0] = (fai0,y)/(fai0,fai0)
    double tmp1=0,tmp2=0;
    for (int i=0;i<m;i++){
        tmp1 += w[i] * y[i];
        tmp2 += w[i];
    }
    a[0] = tmp1 / tmp2;
    //P(x) = a[0] * fai0(x)
    for (int i=0;i<=MAX_n;i++){
        c[i] = a[0] * fai0[i];
    }
    //err = (y,y) - a[0] * (fai0,y)
    double tmp3=0;
    for (int i=0;i<m;i++){
        tmp3 += w[i] * y[i] * y[i];
    }
    err = tmp3 - a[0] * tmp1;

    //STEP 2
    //B = (x*fai0,fai0) / (fai0,fai0) ; fai1(x) = x - B
    double tmp4=0;
    for (int i=0;i<m;i++){
        tmp4 += w[i] * x[i];
    }
    fai1[0] = -(tmp4/tmp2);
    fai1[1] = 1;
    //a[1] = (fai1,y) / (fai1,fai1)
    tmp1 = tmp2 = 0;
    for (int i=0;i<m;i++){
        tmp1 += w[i] * (fai1[0] + x[i]) * y[i];
        tmp2 += w[i] * (fai1[0] + x[i]) * (fai1[0] + x[i]);
    }
    a[1] = tmp1 / tmp2;
    //P(x) += a[1] * fai1(x)
    for (int i=0;i<=MAX_n;i++){
        c[i] += a[1] * fai1[i];
    }
    //err -= a[1] * (fai1,y);
    err -= a[1] * tmp1;

    //STEP 3 4
    int k = 1;
    while ((k<MAX_n)&&(fabs(err)>=(*eps))) {
        k++; //STEP 5
        //STEP 6
        //B = (x*fai1,fai1) / (fai1,fai1)
        //C = (x*fai1,fai0) / (fai0,fai0)
        tmp1 = tmp2 = tmp3 = tmp4 = 0;
        for (int i=0;i<m;i++){
            double tmp5=0,tmp6=0; //tmp5 for fai1, tmp6 for fai0
            for (int j=0;j<=MAX_n;j++){
                tmp5 += fai1[j] * pow(x[i],j);
                tmp6 += fai0[j] * pow(x[i],j);
            }
            tmp1 += x[i] * tmp5 * tmp5;
            tmp2 += tmp5 * tmp5;
            tmp3 += x[i] * tmp5 * tmp6;
            tmp4 += tmp6 * tmp6;
        }
        B = tmp1 / tmp2;
        C = tmp3 / tmp4;
        //fai2(x) = (x-B) * fai1(x) - C * fai0(x)
        fai2[0] = - B * fai1[0] - C * fai0[0];
        for (int i=1;i<=MAX_n;i++){
            fai2[i] = fai1[i-1] - B * fai1[i] - C * fai0[i];
        }
        //a[k] = (fai2,y) / (fai2,fai2)
        tmp1 = tmp2 = tmp3 = tmp4 = 0;
        for (int i=0;i<m;i++){
            double tmp5=0;
            for (int j=0;j<=MAX_n;j++){
                tmp5 += fai2[j] * pow(x[i],j);
            }
            tmp1 += tmp5 * y[i];
            tmp2 += tmp5 * tmp5;
        }
        a[k] = tmp1 / tmp2;
        //P(x) += a[k] * fai2(x)
        for (int i=0;i<=MAX_n;i++){
            c[i] += a[k] * fai2[i];
        }
        //err -= a[k] * (fai2,y)
        err -= a[k] * tmp1;

        //STEP 7
        //fai0(x)=fai1(x),fai1(x)=fai2(x);
        for (int i=0;i<=MAX_n;i++){
            fai0[i] = fai1[i];
            fai1[i] = fai2[i];
        }
    }
    //STEP 8
    *eps = err;
    return k;
}