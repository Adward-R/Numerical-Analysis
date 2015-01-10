//
//  nalab06.cpp
//  nalab06
//
//  Created by Adward on 14/12/22.
//
//

#include <iostream>
#define MAX_N 20
using namespace std;

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn,
                  double a[], double b[], double c[], double d[]);

double S( double t, double Fmax,
         int n, double x[], double a[], double b[], double c[], double d[] );

int main()
{
    int n, Type, m, i;
    double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
    double s0, sn, Fmax, t0, tm, h, t;
    
    while (scanf("%d", &n) != EOF) {
        for (i=0; i<=n; i++)
            scanf("%lf", &x[i]);
        for (i=0; i<=n; i++)
            scanf("%lf", &f[i]);
        scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);
        
        Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
        for (i=1; i<=n; i++)
            printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);
        
        scanf("%lf %lf %d", &t0, &tm, &m);
        h = (tm-t0)/(double)m;
        for (i=0; i<=m; i++) {
            t = t0+h*(double)i;
            printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
        }
        printf("\n");
    }
    
    return 0;
}

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn,double a[], double b[], double c[], double d[]){
    for (int i=0;i<=n;i++) a[i] = f[i];

    double l[n+1],miu[n+1],z[n+1];
    
    if (Type==2) {

        l[0] = 1;
        miu[0] = z[0] = 0;
        for (int i=1;i<n;i++){
            l[i] = 2 * (x[i+1]-x[i-1]) - (x[i]-x[i-1]) * miu[i-1];
            miu[i] = (x[i+1]-x[i]) / l[i];
            z[i] = ( 3 * (f[i+1]-f[i]) / (x[i+1]-x[i]) - 3 * (f[i]-f[i-1]) / (x[i]-x[i-1]) - (x[i]-x[i-1]) * z[i-1] ) / l[i];
        }
        
        l[n] = 1;
        c[n] = z[n] = 0;
        for (int j=n-1;j>=0;j--){
            c[j] = z[j] - miu[j] * c[j+1];
            b[j] = (a[j+1] - a[j]) / (x[j+1] - x[j]) - (x[j+1] - x[j]) * (c[j+1] + 2*c[j]) / 3;
            d[j] = (c[j+1] - c[j]) / (3*(x[j+1] - x[j]));
        }
    }
    else {
        l[0] = 2 * (x[1] - x[0]);
        miu[0] = 0.5;
        z[0] = ( 3 * (a[1] - a[0]) / (x[1] - x[0]) - 3 * s0) / ( 2 * (x[1] - x[0]) );

        for (int i=1;i<n;i++){
            l[i] = 2 * (x[i+1] - x[i-1]) - (x[i] - x[i-1]) * miu[i-1];
            miu[i] = (x[i+1] - x[i]) / l[i];
            z[i] = ( 3 * (f[i+1]-f[i]) / (x[i+1]-x[i]) - 3 * (f[i] - f[i-1]) / (x[i]-x[i-1]) - (x[i]-x[i-1]) * z[i-1] ) / l[i];
        }

        l[n] = (x[n]-x[n-1]) * (2 - miu[n-1]);
        z[n] = (3 * sn - 3 * (a[n]-a[n-1]) / (x[n]-x[n-1]) - (x[n]-x[n-1]) * z[n-1]) / l[n];
        c[n] = z[n];
        
        for (int j=n-1;j>=0;j--){
            c[j] = z[j] - miu[j] * c[j+1];
            b[j] = (a[j+1] - a[j]) / (x[j+1] - x[j]) - (x[j+1]-x[j]) * (c[j+1] + 2 * c[j]) / 3;
            d[j] = (c[j+1] - c[j]) / (3 * (x[j+1]-x[j]));
        }
    }
    for (int i=n;i>=0;i--){
        a[i] = a[i-1];
        b[i] = b[i-1];
        c[i] = c[i-1];
        d[i] = d[i-1];
    }
}

double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] ){
    int j=0;
    while (j<n){
        if (t >= x[j] && t <= x[j+1]) break;
        else j++;
    }
    if (j==n) return Fmax;
    else return a[j+1] + b[j+1] * (t-x[j]) + c[j+1] * (t-x[j]) * (t-x[j]) + d[j+1] * (t-x[j]) * (t-x[j]) * (t-x[j]);
}