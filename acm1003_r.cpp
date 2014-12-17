#include <stdio.h>
#include <math.h>
#include <iostream>

#define Max_size 10000 /* max number of dishes */
void Price( int n, double p[] );
 
int main()
{
    int n, i;
    double p[Max_size];
 
    while (scanf_s("%d", &n)!=EOF) {
       for (i=0; i<n; i++)
           scanf_s("%lf", &p[i]);
       Price(n, p);
       for (i=0; i<n; i++)
           printf("%.2f ", p[i]);
       printf("\n");
    }
 
    return 0;
}

/*void test(double *num, double *p, int n){
int i, j;
for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++)
	   printf("%5.2lf", num[i * n + j]);
	   printf(" ");
	   printf("%5.2lf", p[i]);
	   printf("\n");
}
}

void cal_j(double *num, double *p, int i, int j, int n, double k){
p[j - 1] = p[j - 1] - k * p[i - 1];
for(int m = 0; m < n; m++)
num[(j - 1) * n + m] = num[(j - 1) * n + m] - k * num[(i - 1) * n + m];
}

void Price(int n, double p[]){
double *a;
int i, j, m;
a = new double[n * n];
for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
        if((j - i) == 1 || (i - j) == 1)a[i * n + j] = 0.5;
		else if(i == j)a[i * n + j] = 2;
		     else a[i * n + j] = 0; 
    }
	a[n - 1] = 0.5;
	a[n * (n - 1)] = 0.5;
	//test(a, p, n);
	for(m = 1; m < n - 1; m++){
	   cal_j(a, p, m, m + 1, n, a[m * n + m - 1]/a[(m - 1) * n + m - 1]);
	   //test(a, p, n);
	}
	for(m = 0; m < n - 1; m++)
	   cal_j(a, p, m+1, n, n, a[n * (n - 1) + m]/a[m * n + m]);
	//test(a, p, n);
	   p[n - 1] = p[n - 1]/a[n * n - 1];
	   p[n - 2] = (p[n - 2] - a[(n - 1) * n - 1] * p[n - 1])/a[n * (n - 1) - 2];
	   for(i = n - 3; i >= 0; i--)
       p[i] = (p[i] - a[i * n + n - 1] * p[n - 1] - a[i * n + i + 1] * p[i + 1])/a[i * n + i];  	   
}*/

void cal_j(double *num, double *numl, double *p, int i, int j, int n, double k){
p[j - 1] = p[j - 1] - k * p[i - 1];
if(j < (n - 1)){
   //for(int m = 0; m < 3; m++)
       num[(j - 1 ) * 3] = num[(j - 1) * 3 + 1] - k * num[(i - 1) * 3 + 1];
       num[(j - 1) * 3 + 1] = 0.5;
       num[(j - 1) * 3 + 2] = -k * num[(i - 1) * 3 + 2];
    }
   else {
	   if(j == (n - 1)){
         num[(n - 2) * 3] = 0;
		 num[(n - 2) * 3 + 1] = num[(j - 1) * 3 + 1] - k * num[(i - 1) * 3 + 1];
		 num[(n - 2) * 3 + 2] +=  - k * num[(i - 1) * 3 + 2];
     }
	 else {
		 if(j == n){
	 
	 numl[n - 1] -= k * num[(i - 1) * 3 + 2];
	 numl[i - 1] = 0;
	 if(i != (n - 1))
     numl[i]     -= k * num[(i - 1) * 3 + 1];
	 }
	 }
   }
}

void Price(int n, double p[]){
double *a, *al;
int i, j, m;
a = new double[(n - 1) * 3];
al = new double[n];
/*for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
        if(abs(j - i) == 1)a[i * n + j] = 0.5;
		else if(i == j)a[i * n + j] = 2;
		     else a[i * n + j] = 0; 
    }*/
	a[0] = 2; a[1] = 0.5; a[2] = 0.5;
	for(i = 1; i < n - 1; i++){
	a[i * 3] = 0.5; a[i * 3 + 1] = 2; a[i * 3 + 2] = 0.5;
	}
	for(i = 0; i < n; i++){
	if(i == 0)al[i] = 0.5;
	else if(i == (n - 2))al[i] = 0.5;
	    else if(i == (n - 1)) al[i] = 2;
		      else al[i] = 0;
	}
	
	//test(a, p, n);
	for(m = 1; m < n - 1; m++){
	   cal_j(a, al, p, m, m + 1, n, 0.5/a[3 * (m - 1)]);
	   //test(a, p, n);
	}
	//cal_j(a, al, p, m, m + 1, n, 0.5/a[3 * m + 1]);
	for(m = 0; m < n - 2; m++)
	   cal_j(a, al, p, m+1, n, n, al[m]/a[m * 3]);
	cal_j(a, al, p, m+1, n, n, al[m]/a[m * 3 + 1]);
	//test(a, p, n);
	   p[n - 1] = p[n - 1]/al[n - 1];
	   p[n - 2] = (p[n - 2] - a[(n - 2) * 3 + 2] * p[n - 1])/a[(n - 2) * 3 + 1];
	   for(i = n - 3; i >= 0; i--)
       p[i] = (p[i] - a[i * 3 + 2] * p[n - 1] - a[i * 3 + 1] * p[i + 1])/a[i * 3];  	   
}