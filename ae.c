//ae.c: Cosmic evolution for Einstein-Aether gravity
 
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
//#arbitrary changes 
//Set physical parameters
#define H0 70.0
#define Theta_k 0.0
 
//#define alpha 0.2//1.0
#define kappa 1.0
#define n 1.0
//#define gamma 1.0
 
//Set initial conditions
#define X_init 1.0
#define Y_init 0.8
#define Xdot_init H0
#define Ydot_init 0.0
 
//Set computational parameters
#define max_iter 12000000
#define stepsize -0.00000001 //In FLRW cosmology, we integrate backwards in time
#define filter 10
//#define EPS 0.0001
#define X_infty 1600000
//#define X_max 0.4166666667
 
#define zeta_init 0.0
#define zeta_final 100000000*1.570796327
 
//Define functions 
double Xddot(double X, double Y, double Xdot, double Ydot, double Theta_m, double Theta_Lambda, double alpha, double gamma)
{
    double H1, H2, T1, P, Pp, Ppp;
 
    P=gamma*pow(3.*alpha*Xdot*Xdot/X/X,n);
        Pp=n*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-1);
        Ppp=n*(n-1)*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-2);
 
    H1=2.-alpha*Pp-6.*alpha*alpha*Ppp*Xdot*Xdot/X/X;
    H2=alpha*Pp-6.*alpha*alpha*Ppp*Xdot*Xdot/X/X;
 
    T1 = H2*Xdot*Xdot/(X*X) - P/3. + 2.*Theta_Lambda - Theta_m/(X*X*X);
 
    return X*T1/H1;
}
 
double Yddot(double X, double Y, double Xdot, double Ydot, double Theta_m, double Theta_Lambda, double alpha, double gamma)
{    
    return 1.0;
}
 
double zetadot(double X, double Y)
{
    return -1.0/X;
}
 
double epsilon(double X, double Y, double Xdot, double Ydot, double Theta_m, double Theta_Lambda, double alpha, double gamma)
{
    double G1, P, Pp, Ppp;
 
    P=gamma*pow(3.*alpha*Xdot*Xdot/X/X,n);
        Pp=n*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-1.);
        Ppp=n*(n-1)*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-2.);
    G1=1.-alpha*Pp;
 
    return G1*Xdot*Xdot/(X*X) - Theta_k/(X*X) - Theta_Lambda - Theta_m/(X*X*X)+P/6.;
 
}
 
double Theta_Lambdasolve(double X, double Y, double Xdot, double Ydot, double Theta_m, double alpha, double gamma)
{
    double G1, P, Pp, Ppp;
 
        P=gamma*pow(3.*alpha*Xdot*Xdot/X/X,n);
        Pp=n*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-1.);
        Ppp=n*(n-1)*gamma*pow(3.*alpha*Xdot*Xdot/X/X,n-2.);
        G1=1.-alpha*Pp;
 
    return G1*Xdot*Xdot/(X*X) - Theta_k/(X*X) - Theta_m/(X*X*X)+P/6.;
}
 
double dsolve(double Theta_m, double Theta_Lambda, double alpha, double gamma, double X_max)
{
    double X, Xdot, X1, Xdot1, Y, Ydot, Y1, Ydot1, zeta;
    double tau, tau_check, x, y, z;
    double h, k1, k2, k3, k4, l1, l2, l3, l4, p1, p2, p3, p4, q1, q2, q3, q4, m1, m2, m3, m4;
    int i, i_check;
 
    h=stepsize;
    X=1.0; //X_init;
        Xdot=H0;
        Y=Y_init;
        Ydot=Y_init;
    zeta=zeta_init;
        tau=0;
 
 
    for(i=0; i<max_iter; i++)
    {
        X1=X;
                Y1=Y;
                Xdot1=Xdot;
                Ydot1=Ydot;
                k1=h*Xdot1;
                p1=h*Ydot1;
                l1=h*Xddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
                q1=h*Yddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
        m1=h*zetadot(X1, Y1);
 
                X1=X+k1/2.;
                Y1=Y+p1/2.;
                Xdot1=Xdot+l1/2.;
                Ydot1=Ydot+q1/2.;
                k2=h*Xdot1;
                p2=h*Ydot1;
                l2=h*Xddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
                q2=h*Yddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
        m2=h*zetadot(X1, Y1);
 
                X1=X+k2/2.;
                Y1=Y+p2/2.;
                Xdot1=Xdot+l2/2.;
                Ydot1=Ydot+q2/2.;
                k3=h*Xdot1;
                p3=h*Ydot1;
                l3=h*Xddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
                q3=h*Yddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
        m3=h*zetadot(X1, Y1);
 
                X1=X+k3;
                Y1=Y+p3;
                Xdot1=Xdot+l3;
                Ydot1=Ydot+q3;
                k4=h*Xdot1;
                p4=h*Ydot1;
                l4=h*Xddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
                q4=h*Yddot(X1, Y1, Xdot1, Ydot1, Theta_m, Theta_Lambda, alpha, gamma);
        m4=h*zetadot(X1, Y1);
 
        X = X + (k1 + 2*k2 + 2*k3 + k4)/6.;
                Y = Y + (p1 + 2*p2 + 2*p3 + p4)/6.;
                Xdot = Xdot + (l1 + 2*l2 + 2*l3 + l4)/6.;
                Ydot = Ydot + (q1 + 2*q2 + 2*q3 + q4)/6.;
        zeta = zeta + (m1 + 2*m2 + 2*m3 + m4)/6.;
 
        tau=i*h;
 
        //Break conditions
        if(X>X_infty || X<X_max)
        {
            //printf("%f %f %f %f %f %.10f \n",x,y,z,X,epsilon(X,Y,Xdot,Ydot),zeta);
            break;
        }
        if(i>max_iter)
        {
            break;
        }
 
                /*
                tau_check=(tau/h)/filter;
                if(floor(tau_check)==tau_check)
                {
            printf("%.10f %.10f %.10f %.10f \n",X,epsilon(X,Y,Xdot,Ydot, Theta_m, Theta_Lambda, alpha, gamma),tau, zeta/X);
                }
        */
    }    
    
    return zeta/X;
}

//-----------------------------------------//
//main programme 
int main()
{
    double Theta_msub, Theta_Lambdasub, alphasub, gammasub, dL, mu;
    double Hub = 18.16011608;
	const int SIZE = 580;
	double z_array[SIZE]; // an array to store z values
	double mu_array[SIZE]; // an array to store the dsolved mus
	
	Theta_msub=1470.0;
    alphasub=0.2;
    gammasub=0.0; 
    Theta_Lambdasub=Theta_Lambdasolve(1.0, 1.0, H0, 1.0, Theta_msub, alphasub, gammasub);   

	//read in the data 
    FILE *fptr;
    int i = 0;
   
	fptr = fopen("SCPUnion21z","r");
	
	// check that the file can be opened
	if (!fptr) {
		fprintf(stderr, "Error opening file for reading.\n");
		return 1;
	}

	for (i = 0; i < SIZE; i++) {
		fscanf(fptr, "%lf", &z_array[i]);
	}   

	fclose(fptr);
	
	for (i = 0; i < SIZE; i++) {    
		dL=dsolve(Theta_msub, Theta_Lambdasub, alphasub, gammasub, 1.0/(1+z_array[i]));
		mu_array[i] = 25.0+Hub+5*log10(dL);
	}
    
    //export mu_array into a txt file   
	fptr = fopen("mu.txt","w");
	
	for (i = 0; i < SIZE; i++) {
		fprintf(fptr, "%.10f \n", mu_array[i]);
	}   
   
	fclose(fptr);
	
	
    return 0;
}
