#include<stdio.h>
#include<math.h>
#include<malloc.h>

double dx=0.025,dy=0.025,dt=0.005;
double alpha;
double pi=3.141592653589793238462643383279502884197169399375105820974944592;

int i,j,n;
double *double_star;

double u_exact(double x,double y,double t);
double delta_x2(double **u,int i,int j);
double delta_y2(double **u,int i,int j);
double* chasing(int order,double *d);

int main(){
    alpha=dt/dx/dx;
    int im=(int)(1.0/dx),jm=(int)(1.0/dy),nm=(int)(1.0/dt);
    double ***u,**u_star;
    u=(double ***)malloc((nm+1)*sizeof(double **));
    for(n=0;n<=nm;n++){
        u[n]=(double **)malloc((im+1)*sizeof(double *));
        for(i=0;i<=im;i++){
            u[n][i]=(double *)malloc((jm+1)*sizeof(double));
        }
    }

    u_star=(double **)malloc((im+1)*sizeof(double *));
    for(i=0;i<=im;i++){
        u_star[i]=(double *)malloc((jm+1)*sizeof(double));
    }
    static double x,t,y;
    //设定初始条件
    for(i=1;i<im;i++){
        for(j=1;j<jm;j++){
            y=dy*j;
            x=dx*i;
            t=dt*n;
            u[0][i][j]=20+80*(y-sin(0.5*pi*x)*sin(0.5*pi*y));
        }
    }
    //设定边界条件1
    for(n=0;n<=nm;n++){
        for(j=1;j<jm;j++){
            y=dy*j;
            t=dt*n;
            u[n][0][j]=20+80*y;
            u[n][im][j]=20+80*(y-exp(-0.5*pi*pi*t)*sin(0.5*pi*y));
        }
    }
    //设定边界条件2
    for(n=0;n<=nm;n++){
        for(i=0;i<=im;i++){
            x=dx*i;
            t=dt*n;
            u[n][i][0]=20;
            u[n][i][jm]=20+80*(1-exp(-0.5*pi*pi*t)*sin(0.5*pi*x));
        }
    }
    double d[im];
    for(n=0;n<nm;n++){
        //第一个方向
        for(j=1;j<jm;j++){
            u_star[0][j]=0.5*(u[n][0][j]+u[n+1][0][j])-0.25*alpha*(delta_y2((double **)u[n+1],0,j)-delta_y2((double **)u[n],0,j));
            u_star[im][j]=0.5*(u[n][0][j]+u[n+1][0][j])-0.25*alpha*(delta_y2((double **)u[n+1],im,j)-delta_y2((double **)u[n],im,j));
            //设置d向量
            
            for(i=1;i<im;i++){
                d[i]=(1-alpha)*u[n][i][j]+0.5*alpha*(u[n][i][j+1]+u[n][i][j-1]);
            }
            d[1]+=0.5*alpha*u_star[0][j];
            d[im-1]+=0.5*alpha*u_star[im][j];
            double_star=chasing(im-1,d+1);

            for(i=1;i<im;i++){
                u_star[i][j]=double_star[i-1];
            }
            free(double_star);

        }

        //第二方向
        for(i=1;i<im;i++){
            //设置d
            for(j=1;j<jm;j++){
                d[j]=u_star[i][j]+0.5*alpha*delta_x2((double **)u_star,i,j);
            }
            d[1]+=0.5*alpha*u[n+1][i][0];
            d[im-1]+=0.5*alpha*u[n+1][i][im];
            double_star=chasing(im-1,d+1);
            for(j=1;j<im;j++){
                u[n+1][i][j]=double_star[j-1];
            }
            free(double_star);
        }
    }
    //输出数据
    FILE *fp;
    fp=fopen("dataout.csv","w");
    for(n=0;n<=nm;n++){
        for(i=0;i<=im;i++){
            for(j=0;j<=jm;j++)fprintf(fp,"%g,",u[n][i][j]);
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
    fp=fopen("errorout.csv","w");
    for(n=0;n<=nm;n++){
        for(i=0;i<=im;i++){
            for(j=0;j<=jm;j++)fprintf(fp,"%g,",u[n][i][j]-u_exact(i*dx,j*dy,n*dt));
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
    fp=fopen("exact.csv","w");
    for(n=0;n<=nm;n++){
        for(i=0;i<=im;i++){
            for(j=0;j<=jm;j++)fprintf(fp,"%g,",u_exact(i*dx,j*dy,n*dt));
            fprintf(fp,"\n");
        }
    }
    fclose(fp);

    return 0;
}

double delta_x2(double **u,int i,int j){
    return -u[i][j]*2+u[i-1][j]+u[i+1][j];
}
double delta_y2(double **u,int i,int j){
    return -u[i][j]*2+u[i][j-1]+u[i][j+1];
}
//追赶法,a=1+alpha,b=-0.5*alpha
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

double u_exact(double x,double y,double t){
    return 20+80*(y-exp(-0.5*pi*pi*t)*sin(0.5*pi*x)*sin(0.5*pi*y));
}