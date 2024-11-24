#include<stdio.h>
#include<math.h>
double v(double r);
double f(double r,double E,int l);

int main()
{
int i,l,n,preci=0;
double h=0.003,d=0.1,r[60000],u[60000],w[60000],k[21000],D[6000],NumTD,DenTD,TD,DELTA,Kappa,U[60000],R[60000],r1,r2,U1,U2,R1[2000],R2[2000],E,sigma[21000],pie=3.14,huge,sum[2100],M,dum,dum2,check=0;
double fr,rsqr,Z=1,a1=16.039,a2=2.007,a3=-25.543,a4=4.525,a5=0.961,a6=0.443;
FILE *fp;
FILE *fp2;
FILE *fp3;
FILE *fp4;
fp=fopen("final_bwfn_r_14au.in","w");
fp2=fopen("unormbwfn_total.dat","w");
fp3=fopen("fr_r.dat","w");
fp4=fopen("r_v.dat","w");
//  E is the binding energy of 3p e of Ar
//E=-0.579689295;
E=-0.8; 
M=1;
l=1;
for(n=0;n<=1000;n++)
{
E=E+d;
r[0]=0.001;
r[1]=r[0]+h;
u[0]=pow(r[0],l+1);
u[1]=pow(r[1],l+1);
//huge=fabs(u[1]);

//printf("%f  %f %f \n ",u[i],fabs(u[i]),huge);
w[0]=u[0]*(1-(h*h*f(r[0],E,l))/12);
w[1]=u[1]*(1-(h*h*f(r[1],E,l))/12);
//printf("r[0] is %f,r[1] is %f,u[0] is %f,u[1] is %f",r[0],r[1],u[0],u[1]);

for (i=2;i<=8340;i++)
{
r[i]=r[i-1]+h;
w[i]=-w[i-2]+(2.0*w[i-1])+(h*h*f(r[i-1],E,l)*u[i-1]);
u[i]=w[i]/(1-(h*h*f(r[i],E,l)/12.0));
//if (fabs(u[i]) >huge)
   //huge=fabs(u[i]);
return 0;
}

//printf("r[0] is %lf,r[1] is %lf,r[2]is %lf,u[0] is %lf,u[1] is %lf,u[2] is %lf",r[0],r[1],r[2],u[0],u[1],u[2]);

for(i=0;i<=8340;i++)
{
U[i]=u[i];
R[i]=U[i]/r[i];
//if (n==92)
//{
//fprintf(fp2,"%lf  %lf   %lf  %lf \n" ,E,r[i],U[i],R[i]);
dum=f(r[i],E,l);
//fprintf(fp3,"%lf   %lf \n" ,r[i],dum);
dum2=(Z+a1*exp(-a2*r[i])+a3*r[i]*exp(-a4*r[i])+a5*exp(-a6*r[i]))/r[i];

//fprintf(fp4,"%lf  %lf  \n",r[i],dum2);
//}
}
r1=r[4704];
r2=r[4740];
U1=U[4704];
U2=U[4740];
R1[n]=R[4704];
R2[n]=R[4740];
//if (n==34||n==31||n==32||n==33)
//printf(" %d    %lf    %lf   \n", n, E[n],DELTA*180/4.1415);
 if(n>=1)
	{
	check=R1[n]*R1[n-1];
		if(check<0)
		{
		preci=preci+1;
		E=E-d;
		d=d/10;
		n=0;
		}
	}
	if(preci>8)
	{
	printf("Energy is %lf \n",E);
	printf("Precision is %d",preci);
	break;
	}
 fprintf(fp," %lf  %d  %lf   %lf \n", E,n,R1[n],check);
//fprintf(fp," %lf  %lf  %lf   %lf  %lf  %lf  %lf\n", E,r1,U1,R1,r2,U2,R2);
}
for(i=0;i<=8340;i++)
{
fprintf(fp2,"%lf  %lf   %lf  %lf \n" ,E,r[i],U[i],R[i]);
}
}
// v(r) is confinement potential

double v (double r)
{
double vr;
vr=0; 
return (vr);
}
double f(double r,double E, int l)
{
// Effective atomic no Z,
double fr, vr,rsqr,M=1,Z=1,a1=16.039,a2=2.007,a3=-25.543,a4=4.525,a5=0.961,a6=0.443;
double dum2;
rsqr=r*r;
//Add Coulomb term
dum2=(Z+a1*exp(-a2*r)+a3*r*exp(-a4*r)+a5*exp(-a6*r))/r;
fr=(v(r)-E)*2*M+(l*(l+1)/(rsqr))-(2*M*dum2);
return (fr);
}


