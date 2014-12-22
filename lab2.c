#include <stdio.h>
#include <math.h>
 
#define ZERO 0.000000001 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11   /* Max Polynomial Degree + 1 */
#define N 10000
 
double Polynomial_Root(int n, double c[], double a, double b, double EPS);
int main()
{
    int n;
    double c[MAXN], a, b;
    double EPS = 0.00005;
    int i;
 
    while (scanf("%d", &n)!= EOF){
       for (i=n; i>=0; i--)
           scanf("%lf", &c[i]);
       scanf("%lf %lf", &a, &b);
       printf("%.4lf\n", Polynomial_Root(n, c, a, b, EPS));
    }
 
    return 0;
}


double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    double p,p0=(a+b)/2;
    int i,j;
    double f,f_;
    while (i<=N){
        f=c[0];
        for (j=1;j<=n;j++){
            f+=c[j]*pow(p0,j);
        }
        f_=c[1];
        for (j=2;j<=n;j++){
            f_+=j*c[j]*pow(p0,j-1);
        }
        p=p0-f/f_;
        if (p-p0<EPS&&p-p0>-EPS){
            return p;
        }
        else{
            i++;
            p0=p;
        }
    }
    printf("Failed!\n");
    return 0.0;
}