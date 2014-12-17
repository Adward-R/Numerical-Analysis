#include <stdio.h>
void Series_Sum(double sum[]);
int main(){
    int i;
    double x, sum[3001];
    Series_Sum( sum );
    x = 0.0;
    for (i=0; i<3001; i++)
        printf("%6.2lf %16.12lf\n", x + (double)i * 0.10, sum[i]);
    return 0;
}
/*
void Series_Sum(double sum[]){
    int i,j,n=1222704;
    double x;
    for (i=0;i<3001;i++){
        sum[i]=1;
        x=(double)i * 0.10;
        for (j=1;j<n+1;j++){
            sum[i]+=(1.0-x)/((double)j*((double)j+1.0)*((double)j+x));
        }
    }
}
*/

void Series_Sum(double sum[]){
    int i,j,n=20000;
    double x;
    for (i=0;i<3001;i++){
        x=(double)i * 0.10;
        sum[i]=1.25-0.25*x;
        for (j=1;j<n+1;j++)
            sum[i]+=(2.0-x)*(1.0-x)/((double)j*((double)j+1.0)*((double)j+2.0)*((double)j+x));
    }
}
/*
void Series_Sum(double sum[]){
    int i,j,n=2672;
    double x;
    for (i=0;i<3001;i++){
        x=(double)i * 0.10;
        sum[i]=1.0+0.25*(1.0-x)+(1.0-x)*(2.0-x)/18+(1.0-x)*(2.0-x)*(3.0-x)/96+(1.0-x)*(2.0-x)*(3.0-x)*(4.0-x)/600+(1.0-x)*(2.0-x)*(3.0-x)*(4.0-x)*(5.0-x)/4320;
        for (j=1;j<n+1;j++){
            sum[i]+=(6.0-x)*(5.0-x)*(4.0-x)*(3.0-x)*(2.0-x)*(1.0-x)/(((double)j)*((double)j+1.0)*((double)j+2.0)*((double)j+3.0)*((double)j+4.0)*((double)j+5.0)*((double)j+6.0)*((double)j+x));
        }
    }
    }*/