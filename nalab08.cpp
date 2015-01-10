#include <iostream>
#include <math.h>
using namespace std;
double f0( double x, double l, double t )
{
    return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}
 
double Integral(double a, double b, double (*f)(double x, double y, double z),double eps, double l, double t);
 
int main()
{
    double a=0.0, b, eps=0.005, l, t;
 
    while (scanf("%lf %lf %lf", &l, &b, &t) != EOF)
       printf("%.2f\n", Integral(a, b, f0, eps, l, t));
    return 0;
}


double Romberg(double a, double b, double (*f)(double x, double y, double z),double eps, double l, double t){
    double h = b-a;
    double R[3][100];
    R[1][1] = h/2 * (f(a,l,t)+f(b,l,t));
    int i=2;
    int flag=1;
    do {
        R[2][1] = 0.5 * R[1][1];
        for (int k=1;k<=pow(2,i-2);k++){
            R[2][1] += 0.5 * h * f(a+(k-0.5)*h,l,t);
        }
        for (int j=2;j<=i;j++){
            R[2][j] = R[2][j-1] + (R[2][j-1]-R[1][j-1])/(pow(4,j-1)-1);
        }
        h /= 2;
        if (fabs(R[2][i]-R[1][i-1])<=eps*100) flag=0;
        for (int j=1;j<=i;j++) R[1][j] = R[2][j];
        i++;
    } while (i<=10);
    return R[2][i-1]/100;
}

double Integral(double a, double b, double (*f)(double x, double y, double z),double eps, double l, double t){
    double pi = 3.14159265359;
    double h=b-a;
    double result,re1,re2;
    double w = (2*pi)/t;
    int T = (int) (h/w);

    if (T){
        re1 = Romberg(a,a+w,f0,eps/T,l,t);
        re2 = Romberg(a+w*T,b,f0,eps/T,l,t);
        result = re1*T+re2;
    }
    else{
        result = Romberg(a,b,f0,eps,l,t);
    }
    return result;
}
