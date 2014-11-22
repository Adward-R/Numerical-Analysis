#include <stdio.h>
#define Max_size 10000 /* max number of dishes */
int u,v;

void Price( int n, double p[] );
int main(){
    int n, i;
    double p[Max_size];
    while (scanf("%d", &n)!=EOF) {
       for (i=0; i<n; i++)
           scanf("%lf", &p[i]);
       Price(n, p);
       for (i=0; i<n; i++)
           printf("%.2f ", p[i]);
       printf("\n");
    }
    return 0;
}
void Price(int n, double p[]){
    double a[n][n+1];
    int i,j,k,t;
    double tmp;

    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            a[i][j]=0;
        }
        a[i][i]=4;
        if (i!=0&&i!=n-1){
            a[i][i-1]=a[i][i+1]=1;
        }
    }
    a[0][1]=a[0][n-1]=a[n-1][0]=a[n-1][1]=1;
    for (i=0;i<n;i++){
        a[i][n]=2*p[i];
    }

    for (i=0;i<n-1;i++){
        t=i;
        while (a[t][i]==0 && t<=n-1) p++;
        if (t==n){
            printf("No unique solution exists.\n");
            return;
        }
        else{
            if (t!=i){
                for (j=0;j<=n;j++){
                    tmp=a[t][j];
                    a[t][j]=a[i][j];
                    a[i][j]=tmp;
                }
            }
            for (j=i+1;j<n;j++){
                tmp=a[j][i]/a[i][i];
                for (k=0;k<=n;k++){
                    a[j][k]-=tmp*a[i][k];
                }
            }
        }
    }
    if (a[n-1][n-1]==0){
        printf("No unique solution exists.\n");
        return;
    }

    
   
     printf("\n");
    for (u=0;u<n;u++){
        for (v=0;v<=n;v++){
            printf("%.2f   ",a[u][v]);
        }
        printf("\n");
    }
    
    
    p[n-1]=a[n-1][n]/a[n-1][n-1];

    for (i=n-2;i>=0;i--){
        tmp=0;
        for (j=i+1;j<n;j++){
            tmp+=a[i][j]*p[j];
        }
        p[i]=(a[i][n]-tmp)/a[i][i];
    }
}

