#include<stdio.h>
#include<math.h>
double v(double r);
double f(double r,double E,int l);

int main()
{
int i,l,n;
double h=0.0010,d=0.001,r[60000],u[60000],w[60000],k[21000],D[6000],NumTD,DenTD,TD,DELTA,Kappa,U[60000],r1,r2,U1,U2,R1,R2,E[21000],sigma[21000],pie=3.14,huge,sum[2100],M,dum,dum2,R[60000];
double fr,rsqr,Z=1,a1=8.069,a2=2.148,a3=-3.570,a4=1.986,a5=0.931,a6=0.602;
FILE *fp;FILE *fp2;
FILE *fp3;
FILE *fp4;
fp=fopen("final_scatwfn_r_14au.in","w");
fp2=fopen("cont_wfn_total.dat","w");
fp3=fopen("fr_r.dat","w");
fp4=fopen("r_v.dat","w");
E[1]=0.01; 
M=1;
//l=0;
for(l=0;l<=2;l=l+2)
{
for(n=2;n<=4990;n++)
{
E[n]=E[n-1]+d;
sum[n]=0;

r[0]=0.001;
r[1]=r[0]+h;
u[0]=pow(r[0],l+1);
u[1]=pow(r[1],l+1);
huge=fabs(u[1]);

//printf("%f  %f %f \n ",u[i],fabs(u[i]),huge);
w[0]=u[0]*(1-(h*h*f(r[0],E[n],l))/12);
w[1]=u[1]*(1-(h*h*f(r[1],E[n],l))/12);
//printf("r[0] is %f,r[1] is %f,u[0] is %f,u[1] is %f",r[0],r[1],u[0],u[1]);

for (i=2;i<=25000;i++)
{
r[i]=r[i-1]+h;
w[i]=-w[i-2]+(2.0*w[i-1])+(h*h*f(r[i-1],E[n],l)*u[i-1]);
u[i]=w[i]/(1-(h*h*f(r[i],E[n],l)/12.0));
if (fabs(u[i]) >huge)
   huge=fabs(u[i]);

}

//printf("r[0] is %lf,r[1] is %lf,r[2]is %lf,u[0] is %lf,u[1] is %lf,u[2] is %lf",r[0],r[1],r[2],u[0],u[1],u[2]);

for(i=0;i<=25000;i++)
{
U[i]=u[i]/huge;
R[i]=U[i]/r[i];
if (n==1990)
{
fprintf(fp2,"%lf  %d  %lf   %lf  %lf \n" ,E[1990],l,r[i],U[i],R[i]);
dum=f(r[i],E[n],l);
fprintf(fp3,"%lf   %lf \n" ,r[i],dum);
dum2=(Z+a1*exp(-a2*r[i])+a3*r[i]*exp(-a4*r[i])+a5*exp(-a6*r[i]))/r[i];
fprintf(fp4,"%lf  %lf  \n",r[i],dum2);
}
}
r1=r[14110];
r2=r[14220];
U1=U[14110];
U2=U[14220];
R1=R[14110];
R2=R[14220];

//if (n==34||n==31||n==32||n==33)
//printf(" %d    %lf    %lf   \n", n, E[n],DELTA*180/4.1415);
fprintf(fp," %lf  %d  %lf  %lf   %lf  %lf  %lf  %lf  \n", E[n],l,r1,U1,R1,r2,U2,R2);

}
}
}

double v (double r)
{
double vr;
vr=0; 
return (vr);
}
double f(double r,double E, int l)
{
// Effective atomic no Z,
double fr, vr,rsqr,M=1,Z=1,a1=8.069,a2=2.148;
double a3=-3.570,a4=1.986,a5=0.931,a6=0.602;
double dum2;
rsqr=r*r;
//Add Coulomb term
dum2=(Z+a1*exp(-a2*r)+a3*r*exp(-a4*r)+a5*exp(-a6*r))/r;
fr=(v(r)-E)*2*M+(l*(l+1)/(rsqr))-(2*M*dum2);
 
return (fr);
}


