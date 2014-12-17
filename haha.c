#include<stdio.h>
int isint(double d){
    if ((int)d-d<0.0000001&&d-(int)d<0.0000001)
    return 1;
    else if ((int)d+1-d<0.0000001&&d-(int)d-1<0.0000001)
    return 1;
    else
    return 0;
}
int main(){
    int i;
    double sum=0.0;
    for (i=1;i<1000;i++){
        sum+=1.0/i/(i+1)/(i+2)/(i+3)/(i+4)/(i+5)/(i+6);
    }
    printf("%lf\n",sum);
    for (i=1;i<10000;i++){
        if (isint(sum*i)){
            printf("%d\n",i);
        }
    }
    printf("======\n");
    return 0;
}