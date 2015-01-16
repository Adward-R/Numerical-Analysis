#include <stdio.h>
#include <math.h>
 
#define ZERO 0.000000001 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11   /* Max Polynomial Degree + 1 */
 
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

double cntPolynomial(int m, double c[], double x){
    double p = c[m];
    for (int i=m; i>0; i--) p = p*x + c[i-1];
    return p;
}


int newton(double x0, double *r, double c[], double d[], int n, double a, double b, double eps){
    double cntPolynomial(int m, double c[], double x);
    int i = 1; 
    double p0, p, fp, fq ;

    p = p0 = x0;
    while (i < 1000) {
        fp = cntPolynomial(n, c, p0);
        fq = cntPolynomial(n-1, d, p0);
        if ((fabs(fq)<ZERO) && fabs(fp)>ZERO) return 0;
        if(fabs(fq) == 0){
            *r = p;
            return 1;
        }
        p = ((p0 * fq) - fp)/fq;
        p0 = p;
        i++;
    }
    if(p<a || p>b) return 0;
    *r = p;
    return 1;
}

double Polynomial_Root(int n,double c[], double a, double b, double eps){
    int newton(double x0, double *r,double c[], double d[], int n, double a, double b, double eps);
    double d[MAXN];
    int i;
    double root = 0.0;

    for (i=n-1; i>=0; i--) {
        d[i] = (i+1)*c[i+1];
    }
    if (a>b) {
        root = a; a = b; b = root;
    }
    double average = 0.0;
    int k = 0;
    int t = newton(a, &root, c, d, n, a, b, eps);
    if (t == 1){
        average += root;
        k++;
    }
    int s = newton(b, &root, c, d, n, a, b, eps);
    if (s == 1){
        average += root;
        k++;
    }
    int u = newton((a+b)*0.5, &root, c, d, n, a, b, eps);
    if (u == 1){
        average += root;
        k++;
    }
    average /= k;

    if (fabs(root)<eps) return fabs(root);
    else return root;
}