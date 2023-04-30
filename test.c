#include<stdio.h>
#include<math.h>
#include<malloc.h>

extern double* chasing(int order,double *d);
double alpha=1;
int main(){
    double u[4][4]={{2,1,0,0},{1,2,1,0},{0,1,2,1},{0,0,1,2}};
    double d[4]={1,2,3,4};
    double *x;
    
    x=chasing(4,d);
    return 0;

}

double* chasing(int order,double *d){
        double a=1+alpha,b=-0.5*alpha;
        double *x;
        x=(double *)malloc(order*sizeof(double));
        double uk[order],vk[order],yk[order];
        int i;
        //追
        uk[0]=a;
        vk[0]=b/uk[0];
        yk[0]=d[0]/uk[0];
        for(i=1;i<order;i++){
            uk[i]=a-b*vk[i-1];
            vk[i]=b/uk[i];
            yk[i]=(d[i]-b*yk[i-1])/uk[i];
        }
        //赶
        x[order-1]=yk[order-1];
        for(i=order-2;i>=0;i--){
            x[i]=yk[i]-vk[i]*x[i+1];
        }
        return x;
}