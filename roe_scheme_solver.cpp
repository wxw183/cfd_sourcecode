#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>

/*#include<iostream>
using namespace std;*/

double exact_u(double x);
double exact_p(double x);
double exact_rho(double x);

double max(double a,double b,double c){
    if(a>=b&&a>=c)return(a);
    else if(b>=c)return(b);
    else return(c);
}

int main(){
    int i,j;
    FILE *fp;
    fp=fopen("velocity.csv","w");
    fclose(fp);
    fp=fopen("pressure.csv","w");
    fclose(fp);
    fp=fopen("density.csv","w");
    fclose(fp);
    //探针
    fp=fopen("check.csv","w");
    fclose(fp);
    int net=200;//net为单元数，net+1为单元边界数
    int n=0;//时间步
    double cfl=0.4;//cfl数
    double dt=1.0/net*cfl;
    double c_v=717.5;
    double R=287.0;
    double T_0=300.0;
    double gamma=1.4;
    double *rho,*p,*u,*u_next,*h,*e,*c,*x,*T;//rho为密度，p为压强，u为速度,加next为下一时间步速度，h为焓，e为能量，c为波速，x为横坐标，均为单元上的值
    double *rho_c,*p_c,*u_c,*h_c,*e_c,*c_c,*x_c,*T_c;//为单元边界上的值，rho_c[i]=rho(i-1/2)
    //定义单元上的值
    rho=(double *)malloc(net*sizeof(double));
    p=(double *)malloc(net*sizeof(double));
    u=(double *)malloc(net*sizeof(double));
    u_next=(double *)malloc(net*sizeof(double));
    h=(double *)malloc(net*sizeof(double));
    e=(double *)malloc(net*sizeof(double));
    T=(double *)malloc(net*sizeof(double));
    c=(double *)malloc(net*sizeof(double));
    x=(double *)malloc(net*sizeof(double));

    //定义单元边界上的特征值和特征向量
    double *deltaw_1,*deltaw_2,*deltaw_3;
    deltaw_1=(double *)malloc((net+1)*sizeof(double));
    deltaw_2=(double *)malloc((net+1)*sizeof(double));
    deltaw_3=(double *)malloc((net+1)*sizeof(double));

    //定义单元边界上的值
    rho_c=(double *)malloc((net+1)*sizeof(double));
    p_c=(double *)malloc((net+1)*sizeof(double));
    u_c=(double *)malloc((net+1)*sizeof(double));
    h_c=(double *)malloc((net+1)*sizeof(double));
    e_c=(double *)malloc((net+1)*sizeof(double));
    T_c=(double *)malloc((net+1)*sizeof(double));
    c_c=(double *)malloc((net+1)*sizeof(double));
    x_c=(double *)malloc((net+1)*sizeof(double));

    //定义单元边界上的特征值和特征向量
    double *R_1[3],*R_2[3],*R_3[3],*lambda_1,*lambda_2,*lambda_3,*lambda_1c,*lambda_2c,*lambda_3c;
    for(int i=0;i<3;i++){
        R_1[i]=(double *)malloc((net+1)*sizeof(double));
        R_2[i]=(double *)malloc((net+1)*sizeof(double));
        R_3[i]=(double *)malloc((net+1)*sizeof(double));
    }
    lambda_1c=(double *)malloc((net+1)*sizeof(double));
    lambda_2c=(double *)malloc((net+1)*sizeof(double));
    lambda_3c=(double *)malloc((net+1)*sizeof(double));
    lambda_1=(double *)malloc(net*sizeof(double));
    lambda_2=(double *)malloc(net*sizeof(double));
    lambda_3=(double *)malloc(net*sizeof(double));

    double *epsilon;
    epsilon=(double *)malloc((net+1)*sizeof(double));

    //定义列向量U和F
    double  *U[3],*F[3],*U_c[3],*F_c[3];
    for(int i=0;i<3;i++){
        U[i]=(double *)malloc(net*sizeof(double));
        F[i]=(double *)malloc(net*sizeof(double));
        U_c[i]=(double *)malloc((net+1)*sizeof(double));
        F_c[i]=(double *)malloc((net+1)*sizeof(double));
    }

    for(int i=0;i<net;i++){
        x[i]=1.0/net*(i-0.5);
    }
    //赋初始条件
    for(int i=0;i<net;i++){
        if(x[i]<=0.5){
            rho[i]=1.0;
            p[i]=1.0;
            u[i]=0.0;
        }
        else{
            rho[i]=0.125;
            p[i]=0.1;
            u[i]=0.0;
        }
    }
    
    step1:
    //计算e和h
    for(int i=0;i<net;i++){
        e[i]=0.5*u[i]*u[i]+c_v*p[i]/rho[i]/R;
        h[i]=e[i]+p[i]/rho[i];
    }
    //计算列向量U,F
    for(int i=0;i<net;i++){
        U[0][i]=rho[i];
        U[1][i]=rho[i]*u[i];
        U[2][i]=rho[i]*e[i];
        F[0][i]=rho[i]*u[i];
        F[1][i]=rho[i]*u[i]*u[i]+p[i];
        F[2][i]=rho[i]*u[i]*h[i];
    }
    
    //计算单元边界上的平均值
    for(int i=1;i<net;i++){
        rho_c[i]=sqrt(rho[i-1]*rho[i]);
        u_c[i]=(sqrt(rho[i-1])*u[i-1]+sqrt(rho[i])*u[i])/(sqrt(rho[i-1])+sqrt(rho[i]));
        h_c[i]=(sqrt(rho[i-1])*h[i-1]+sqrt(rho[i])*h[i])/(sqrt(rho[i-1])+sqrt(rho[i]));
        c_c[i]=sqrt((gamma-1)*(h_c[i]-0.5*u_c[i]*u_c[i]));
    }
    //施加边界条件
    rho_c[0]=rho[0];
    u_c[0]=0.0;
    h_c[0]=h[0];
    c_c[0]=sqrt((gamma-1)*(h_c[0]-0.5*u_c[0]*u_c[0]));
    rho_c[net]=rho[net-1];
    u_c[net]=0.0;
    h_c[net]=h[net-1];
    c_c[net]=sqrt((gamma-1)*(h_c[net]-0.5*u_c[net]*u_c[net]));
    
    //施加边界条件(这个是错误的)
    /*u_c[0]=-u_c[1];
    u_c[net]=-u_c[net-1];//边界内移
    rho_c[0]=rho_c[1];
    rho_c[net]=rho_c[net-1];
    c_c[0]=c_c[1];
    c_c[net]=c_c[net-1];
    h_c[0]=h_c[1];
    h_c[net]=h_c[net-1];*/
    //定义单元边界上的特征值和特征向量
    

    for(int i=0;i<net;i++){
        lambda_1[i]=u[i];
        lambda_2[i]=u[i]+c[i];
        lambda_3[i]=u[i]-c[i];
    }
    for(int i=0;i<net+1;i++){
        lambda_1c[i]=u_c[i];
        lambda_2c[i]=u_c[i]+c_c[i];
        lambda_3c[i]=u_c[i]-c_c[i];
        R_1[0][i]=1.0;
        R_1[1][i]=u_c[i];
        R_1[2][i]=0.5*u_c[i]*u_c[i];
        R_2[0][i]=0.5*rho_c[i]/c_c[i];
        R_2[1][i]=0.5*rho_c[i]/c_c[i]*(u_c[i]+c_c[i]);
        R_2[2][i]=0.5*rho_c[i]/c_c[i]*(h_c[i]+u_c[i]*c_c[i]);
        R_3[0][i]=-0.5*rho_c[i]/c_c[i];
        R_3[1][i]=-0.5*rho_c[i]/c_c[i]*(u_c[i]-c_c[i]);
        R_3[2][i]=-0.5*rho_c[i]/c_c[i]*(h_c[i]-u_c[i]*c_c[i]);
    }
    

    //计算波的幅值
    
    for(int i=1;i<net;i++){
        deltaw_1[i]=rho[i]-rho[i-1]-(p[i]-p[i-1])/c_c[i]/c_c[i];
        deltaw_2[i]=u[i]-u[i-1]+(p[i]-p[i-1])/c_c[i]/rho_c[i];
        deltaw_3[i]=u[i]-u[i-1]-(p[i]-p[i-1])/c_c[i]/rho_c[i];

    }
    deltaw_1[0]=0;
    deltaw_2[0]=2*u[0];
    deltaw_3[0]=2*u[0];
    deltaw_1[net]=0;
    deltaw_2[net]=-2*u[net-1];
    deltaw_3[net]=-2*u[net-1];

    
    //熵修正
    
    for(int i=1;i<net;i++){
        epsilon[i]=max(0,lambda_1c[i]-lambda_1[i-1],lambda_1c[i]-lambda_1[i]);
        if(fabs(lambda_1c[i])>=epsilon[i]){
            lambda_1c[i]=fabs(lambda_1c[i]);
        }
        else{
            lambda_1c[i]=0.5*(lambda_1c[i]*lambda_1c[i]/epsilon[i]+epsilon[i]);
        }
    }
    epsilon[0]=max(0,lambda_1c[0]+lambda_1[0],lambda_1c[0]-lambda_1[0]);
    if(fabs(lambda_1c[0])>=epsilon[0]){
            lambda_1c[0]=fabs(lambda_1c[0]);
        }
        else{
            lambda_1c[0]=0.5*(lambda_1c[0]*lambda_1c[0]/epsilon[0]+epsilon[0]);
        }
    epsilon[net]=max(0,lambda_1c[net]-lambda_1[net-1],lambda_1c[net]+lambda_1[net-1]);
    if(fabs(lambda_1c[net])>=epsilon[net]){
            lambda_1c[net]=fabs(lambda_1c[net]);
        }
        else{
            lambda_1c[net]=0.5*(lambda_1c[net]*lambda_1c[net]/epsilon[net]+epsilon[net]);
        }
    
    for(int i=1;i<net;i++){
        epsilon[i]=max(0,lambda_2c[i]-lambda_2[i-1],lambda_2c[i]-lambda_2[i]);
        if(fabs(lambda_2c[i])>=epsilon[i]){
            lambda_2c[i]=fabs(lambda_2c[i]);
        }
        else{
            lambda_2c[i]=0.5*(lambda_2c[i]*lambda_2c[i]/epsilon[i]+epsilon[i]);
        }
    }
    epsilon[0]=max(0,lambda_2c[0]+u[0]-c[0],lambda_2c[0]-lambda_2[0]);
    if(fabs(lambda_2c[0])>=epsilon[0]){
            lambda_2c[0]=fabs(lambda_2c[0]);
        }
        else{
            lambda_2c[0]=0.5*(lambda_2c[0]*lambda_2c[0]/epsilon[0]+epsilon[0]);
        }
    epsilon[net]=max(0,lambda_1c[net]-lambda_1[net-1],lambda_1c[net]+u[net-1]-c[net-1]);
    if(fabs(lambda_2c[net])>=epsilon[net]){
            lambda_2c[net]=fabs(lambda_2c[net]);
        }
        else{
            lambda_2c[net]=0.5*(lambda_2c[net]*lambda_2c[net]/epsilon[net]+epsilon[net]);
        }

    for(int i=1;i<net;i++){
        epsilon[i]=max(0,lambda_3c[i]-lambda_3[i-1],lambda_3c[i]-lambda_3[i]);
        if(fabs(lambda_3c[i])>=epsilon[i]){
            lambda_3c[i]=fabs(lambda_3c[i]);
        }
        else{
            lambda_3c[i]=0.5*(lambda_3c[i]*lambda_3c[i]/epsilon[i]+epsilon[i]);
        }
    }
    epsilon[0]=max(0,lambda_3c[0]+u[0]+c[0],lambda_3c[0]-lambda_3[0]);
    if(fabs(lambda_3c[0])>=epsilon[0]){
            lambda_3c[0]=fabs(lambda_3c[0]);
        }
        else{
            lambda_3c[0]=0.5*(lambda_3c[0]*lambda_3c[0]/epsilon[0]+epsilon[0]);
        }
    epsilon[net]=max(0,lambda_3c[net]-lambda_3[net-1],lambda_3c[net]+u[net-1]+c[net-1]);
    if(fabs(lambda_3c[net])>=epsilon[net]){
            lambda_3c[net]=fabs(lambda_3c[net]);
        }
        else{
            lambda_3c[net]=0.5*(lambda_3c[net]*lambda_3c[net]/epsilon[net]+epsilon[net]);
        }
    
    //计算数值通量
    for(int i=1;i<net;i++){
        for(j=0;j<3;j++){
            F_c[j][i]=0.5*(F[j][i-1]+F[j][i])-0.5*(fabs(lambda_1c[i])*deltaw_1[i]*R_1[j][i])-0.5*(fabs(lambda_2c[i])*deltaw_2[i]*R_2[j][i])-0.5*(fabs(lambda_3c[i])*deltaw_3[i]*R_3[j][i]);
        }

    }
    i=0;
    for(j=0;j<3;j++){
            F_c[j][i]=-0.5*(fabs(lambda_1c[i])*deltaw_1[i]*R_1[j][i])-0.5*(fabs(lambda_2c[i])*deltaw_2[i]*R_2[j][i])-0.5*(fabs(lambda_3c[i])*deltaw_3[i]*R_3[j][i]);
        }
    j=1;
    F_c[j][i]=0.5*(2*F[j][i])-0.5*(fabs(lambda_1c[i])*deltaw_1[i]*R_1[j][i])-0.5*(fabs(lambda_2c[i])*deltaw_2[i]*R_2[j][i])-0.5*(fabs(lambda_3c[i])*deltaw_3[i]*R_3[j][i]);
    i=net;
    for(j=0;j<3;j++){
            F_c[j][i]=-0.5*(fabs(lambda_1c[i])*deltaw_1[i]*R_1[j][i])-0.5*(fabs(lambda_2c[i])*deltaw_2[i]*R_2[j][i])-0.5*(fabs(lambda_3c[i])*deltaw_3[i]*R_3[j][i]);
        }
    j=1;
    F_c[j][i]=0.5*(2*F[j][i-1])-0.5*(fabs(lambda_1c[i])*deltaw_1[i]*R_1[j][i])-0.5*(fabs(lambda_2c[i])*deltaw_2[i]*R_2[j][i])-0.5*(fabs(lambda_3c[i])*deltaw_3[i]*R_3[j][i]);
    /*F_c[0][0]=rho_c[0]*u_c[0];
    F_c[1][0]=rho_c[0]*u_c[0]*u_c[0]+p_c[0];
    F_c[2][0]=rho_c[0]*u_c[0]*h_c[0];
    F_c[0][net+1]=rho_c[net+1]*u_c[net+1];
    F_c[1][net+1]=rho_c[net+1]*u_c[net+1]*u_c[net+1]+p_c[net+1];
    F_c[2][net+1]=rho_c[net+1]*u_c[net+1]*h_c[net+1];*/

    //探针
    fp=fopen("check.csv","a");
    for(int i=0;i<net+1;i++){
        fprintf(fp,"%lg,",rho_c[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
    
    //计算下一时间步的U
    for(int i=0;i<net;i++){
        for(j=0;j<3;j++){
            U[j][i]+=-dt*net*(F_c[j][i+1]-F_c[j][i]);
        }

    }
    
    //反算出下一时间步温度压力密度速度并输出
    for(int i=0;i<net;i++){
        rho[i]=U[0][i];
        u[i]=U[1][i]/rho[i];
        e[i]=U[2][i]/rho[i];
        T[i]=(e[i]-0.5*u[i]*u[i])/c_v;
        p[i]=rho[i]*T[i]*R;
    }
    fp=fopen("density.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",rho[i]);
    }
    fprintf(fp,"%lg\n",rho[net-1]);
    fclose(fp);
    fp=fopen("velocity.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",u[i]);
    }
    fprintf(fp,"%lg\n",u[net-1]);
    fclose(fp);
    fp=fopen("pressure.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",p[i]);
    }
    fprintf(fp,"%lg\n",p[net-1]);
    fclose(fp);
    n++;
    while(n<0.25/dt)goto step1;
    //在文件后面输出0.25秒时解析解
    fp=fopen("pressure.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",exact_p((i+0.5)/net));
    }
    i=net-1;
    fprintf(fp,"%lg\n",exact_p((i+0.5)/net));
    fclose(fp);
    fp=fopen("density.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",exact_rho((i+0.5)/net));
    }
    i=net-1;
    fprintf(fp,"%lg\n",exact_rho((i+0.5)/net));
    fclose(fp);
    fp=fopen("velocity.csv","a");
    for(int i=0;i<net-1;i++){
        fprintf(fp,"%lg,",exact_u((i+0.5)/net));
    }
    i=net-1;
    fprintf(fp,"%lg\n",exact_u((i+0.5)/net));
    fclose(fp);
    //计算密度误差
    double sum=0.0;
    for(int i=0;i<net;i++){
        sum+=pow(exact_rho((i+0.5)/net)-rho[i],2)/net;
    }
    fp=fopen("out.txt","a");
    fprintf(fp,"net=%d,error=%lg\n",net,sum);
    fclose(fp);


    return 0;
}

double exact_p(double x){
    if(x<0.20419601084502)return(1.0);
    else if(x<=0.482431796859703)return(pow(1.11505141824283-0.563436169819009*x,7));
    else if(x<=0.938038933007542)return(0.303130178050645);
    else return(0.1);
}

double exact_rho(double x){
    if(x<0.20419601084502)return(1.0);
    else if(x<=0.482431796859703)return(pow(1.11505141824283-0.563436169819009*x,5));
    else if(x<=0.731863155012236)return(0.426319428178492);
    else if(x<=0.938038933007542)return(0.265573711705307);
    else return(0.125);
}

double exact_u(double x){
    if(x<=0.20419601084502)return(0.0);
    else if(x<=0.482431796859703)return(-0.680653369483392+10.0*x/3.0);
    else if(x<=0.938038933007542)return(0.927452620048944);
    else return(0.0);
}


