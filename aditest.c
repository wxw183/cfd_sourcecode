#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
#include<stdlib.h>
#include<omp.h>

double dx,dy,dt;
double alpha;
double pi=3.141592653589793238462643383279502884197169399375105820974944592;

int i,j,n;


double u_exact(double x,double y,double t);
double delta_x2(double **u,int i,int j);
double delta_y2(double **u,int i,int j);
double* chasing(int order,double *d);
double error_norms(double ***u,int im,int jm);
double epsilon_1(double **u_num,int im,int jm);

int main(int argc,char *argv[]){
    
    int core=atof(argv[4]);
    dx=dy=atof(argv[1]);
    dt=atof(argv[2]);
    double phisical_time=atof(argv[3]);
    //int core=atof(argv[4]);
    if(argc!=5){
        printf("参数数量错误！\n");
        exit(1);
    }
    
    FILE *fp;
    fp=fopen("exact.csv","w");
    fclose(fp);
    fp=fopen("dataout.csv","w");
    fclose(fp);
    fp=fopen("errorout.csv","w");
    fclose(fp);
    fp=fopen("eh.csv","w");
    fclose(fp);
    static double x,t,y;
    alpha=dt/dx/dx;
    int im=(int)(1.0/dx),jm=(int)(1.0/dy),nm=(int)(phisical_time/dt),count;
    double ***u,**u_star;
    u=(double ***)malloc((2)*sizeof(double **));
    
    u[0]=(double **)malloc((im+1)*sizeof(double *));
    for(i=0;i<=im;i++){
        u[0][i]=(double *)malloc((jm+1)*sizeof(double));
    }
    

    u_star=(double **)malloc((im+1)*sizeof(double *));
    for(i=0;i<=im;i++){
        u_star[i]=(double *)malloc((jm+1)*sizeof(double));
    }
    
 //设定初始条件
    for(i=1;i<im;i++){
        for(j=1;j<jm;j++){
            y=dy*j;
            x=dx*i;
            t=dt*n;
            u[0][i][j]=20+80*(y-sin(0.5*pi*x)*sin(0.5*pi*y));
        }
    }

    for(n=0;n<=nm;n++){

        u[1]=(double **)malloc((im+1)*sizeof(double *));
        for(i=0;i<=im;i++){
            u[1][i]=(double *)malloc((jm+1)*sizeof(double));
        }
           
        //设定边界条件1
        for(count=0;count<2;count++){
            for(j=1;j<jm;j++){
                y=dy*j;
                t=dt*(count+n);
                u[count][0][j]=20+80*y;
                u[count][im][j]=20+80*(y-exp(-0.5*pi*pi*t)*sin(0.5*pi*y));
            }
        }
        //设定边界条件2
        for(count=0;count<2;count++){
            for(i=0;i<=im;i++){
                x=dx*i;
                t=dt*(n+count);
                u[count][i][0]=20;
                u[count][i][jm]=20+80*(1-exp(-0.5*pi*pi*t)*sin(0.5*pi*x));
            }
        }
        
        
        //第一个方向
        #pragma omp parallel for private(j) num_threads(core)
        for(j=1;j<jm;j++){
            u_star[0][j]=0.5*(u[0][0][j]+u[1][0][j])-0.25*alpha*(delta_y2((double **)u[1],0,j)-delta_y2((double **)u[0],0,j)); 
            u_star[im][j]=0.5*(u[0][im][j]+u[1][im][j])-0.25*alpha*(delta_y2((double **)u[1],im,j)-delta_y2((double **)u[0],im,j));
            //设置d向量
            double d[im];
            
            for(int i=1;i<im;i++){
                d[i]=(1-alpha)*u[0][i][j]+0.5*alpha*(u[0][i][j+1]+u[0][i][j-1]);
            }
            
            d[1]+=0.5*alpha*u_star[0][j];
            d[im-1]+=0.5*alpha*u_star[im][j];
            double *double_star;
            double_star=chasing(im-1,d+1);
            
            for(int i=1;i<im;i++){
                u_star[i][j]=double_star[i-1];
            }
            //printf("%d\t",j);
            free(double_star);  
        }
        
        //第二方向
        #pragma omp parallel for private(i) num_threads(core)
        for(i=1;i<im;i++){
            //设置d
            double d[im];
            
            for(int j=1;j<jm;j++){
                d[j]=u_star[i][j]+0.5*alpha*delta_x2((double **)u_star,i,j);
            }
            d[1]+=0.5*alpha*u[1][i][0];
            d[im-1]+=0.5*alpha*u[1][i][im];
            double *double_star;
            double_star=chasing(jm-1,d+1);
            
            for(int j=1;j<jm;j++){
                u[1][i][j]=double_star[j-1];
            }
            free(double_star);
        }
        
    //输出数据
    
       
    
       
    
        
        

        
        fprintf(stdout,"%g,%g\n",dt*n,log10(error_norms(u,im,jm)));
        

        

        for(i=0;i<=im;i++){
            free(u[0][i]);
        }
        free(u[0]);
        u[0]=u[1];
    }


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

double error_norms(double ***u,int im,int jm){
    double emax=0.0;
    for(i=0;i<im;i++){
        for(j=0;j<jm;j++){
            if(emax<fabs(u[1][i][j]-u[0][i][j]))emax=fabs(u[1][i][j]-u[0][i][j]);
        }
    }
    return emax;
}

double epsilon_1(double **u_num,int im,int jm){
    double emax=0.0;
    for(i=0;i<im;i++){
        for(j=0;j<jm;j++){
            if(emax<fabs(u_num[i][j]-u_exact(i*dx,j*dy,n*dt)))emax=fabs(u_num[i][j]-u_exact(i*dx,j*dy,n*dt));
        }
    }
    return emax;
}