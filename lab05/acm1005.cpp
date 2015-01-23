#include <stdio.h>
#include <iostream>
#include <iomanip>
#define MAX_SIZE 100

using namespace std; 


int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[],
double TOL, int MAXN);
 
int main()
{
    int n, MAXN, m, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
    double lambda, TOL;
 
    while (scanf_s("%d", &n) != EOF) {
       for (i=0; i<n; i++)
           for (j=0; j<n; j++)
              scanf_s("%lf", &a[i][j]);
       scanf_s("%lf %d", &TOL, &MAXN);
       scanf_s("%d", &m);
       for (i=0; i<m; i++) {
           scanf_s("%lf", &lambda);
           for (j=0; j<n; j++)
              scanf_s("%lf", &v[j]);
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

double Fabs(double a){
    if (a >= 0) {
        return a;
    }
	else {
        return -a;
    }
}

//calculate [A-qI] during the process of LU factorization may reduce the error
int factorization(int n, double a[][MAX_SIZE], double **L, double **U, double q, const double ZERO, double TOL){
    double tem = 0.0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			if(i == j)
		        L[i][j] = U[i][j] = 1.0;
			else
				L[i][j] = U[i][j] = 0;
	}

	U[0][0] = a[0][0] - q;
	if(U[0][0] == 0)return 0;

	for(int j = 1; j < n; j++){
	    U[0][j] = a[0][j];
		L[j][0] = a[j][0] / U[0][0];
	}
	for(int i = 1; i < (n - 1); i++){
	    tem = 0.0;
	    for(int k = 0; k < i; k++)
		    tem += L[i][k] * U[k][i];
		tem = a[i][i] - tem - q;
		if(Fabs(tem) < TOL)tem = 0;
		U[i][i] = tem;
		if(U[i][i] == 0)return 0;

		for(int j = i + 1; j < n; j++){
		    tem = 0.0;
	        for(int k = 0; k < i; k++)
		        tem += L[i][k] * U[k][j];
            tem = a[i][j] - tem;
		    U[i][j] = tem;
		    tem = 0.0;
	        for(int k = 0; k < i; k++)
		        tem += L[j][k] * U[k][i];
		    tem = a[j][i] - tem;
		    L[j][i] = tem / U[i][i];
		}
	}
	tem = 0.0;
	for(int k = 0; k < (n - 1); k++)
		tem += L[n-1][k] * U[k][n-1];
	tem = a[n-1][n-1] - tem - q;
	if(Fabs(tem) < TOL)tem = 0;
	U[n-1][n-1] = tem;
	if(U[n-1][n-1] == 0)return 0;

	if(Fabs(L[n-1][n-1]) < TOL)return 0;

    return 1;
}

int cal_after_LU(int n, double a[][MAX_SIZE], double **L, double **U, double v[], double q, const double ZERO, int first, double TOL){
    int re = -1; 
	double tem = 0.0;
	if (first == 1){
	    re = factorization(n, a, L, U, q, ZERO, TOL);
	    if (re == 0){
            return 0;
        }
	}
	v[0] = v[0];
	for(int i = 1; i < n; i++){
	    tem = 0.0;
		for(int j = 0; j < i; j++)
		    tem += L[i][j] * v[j];
		tem = v[i] - tem;
	    v[i] = tem;	
	}
	v[n-1] = v[n-1] / U[n-1][n-1];
	for(int i = n - 2; i >= 0; i--){
	    tem = 0.0;
		for(int j = i + 1; j < n; j++)
		    tem += U[i][j] * v[j];
		tem = v[i] - tem;
	    v[i] = tem / U[i][i];	
	}

	return 1;
}


/**/
int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN){
    const double ZERO = 0.000000001;
    int k = 1;
	int reLU;
	double *l;
	double err;
	double tem;
	double max;
	double **L, **U;
	reLU = -1;
	L = new double*[n];
	U = new double*[n];
	for(int i = 0; i < n; i++){
	    L[i] = new double[n];
		U[i] = new double[n];
	}
	l = new double[n];
	double q = *lambda;
	max = Fabs(v[0]);
	int index = 0;
	double u = 0.0;
	for(int i = 1; i < n; i++){
	    if(Fabs(v[i]) > max){
            max = Fabs(v[i]);
        	index = i;
        }			
	}
	tem = v[index];
	for(int i = 0; i < n; i++)
	    v[i] /= tem;
	while(k <= MAXN){
		for(int i = 0; i < n; i++){
		    l[i] = v[i];
		}
		//get new v
		if(reLU == -1){
		    reLU = cal_after_LU(n, a, L, U, v, q, ZERO, 1, TOL);
			if(reLU == 0){
				*lambda = q;
				return -1;
			}
        }
		else if(reLU == 1){
			cal_after_LU(n, a, L, U, v, q, ZERO, 0, TOL);
		}

		u = v[index];
		max = u;
		for(int i = 0; i < n; i++){
			if(Fabs(v[i]) > Fabs(max)){
				index = i;
				max = v[i];
			}
		}

		for(int i = 0; i < n; i++){
			v[i] = v[i] / max;
		}
        err = 0;
		for(int i = 0; i < n; i++){
		    if(Fabs(l[i]-v[i]) > err){
			    err = Fabs(l[i]-v[i]);
			}
		}
		
		if(err < TOL){
			u = 1.0 / u + q;
			*lambda = u;
			if(Fabs(*lambda) < ZERO)*lambda = 0;
			for(int i = 0; i < n; i++)
				if(Fabs(v[i]) < ZERO)v[i] = 0;
			return 1;
		}
		k++;
	}
	return 0;
}