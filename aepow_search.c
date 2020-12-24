//aepow.c: Cosmic evolution for Einstein-Aether gravity, power law for the kinetic term

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//Set physical parameters
#define H0 70.0
#define Theta_k 0.0

//#define alpha 0.2//1.0
#define kappa 1.0
#define n 1.0
//#define gamma_top -0.1*5.0*H0*H0
#define gamma_top 0.0
#define Theta_Lambda_top 0.7*H0*H0

//Set initial conditions
#define X_init 1.0
#define Y_init 0.8
#define Xdot_init H0
#define Ydot_init 0.0

//Set computational parameters
#define max_iter 120000000
#define stepsize -0.000000001 //In FLRW cosmology, we integrate backwards in time
#define filter 1000
#define X_infty 2.0
#define EPS 0.01

#define zeta_init 0.0
#define zeta_final 100000000*1.570796327


double Xddot(double X, double Y, double Xdot, double Ydot, double Theta_m, double Theta_Lambda, double alpha, double gamma)
{
	double H, T1;

	H=Xdot/X;

	T1=2.+n*(1-2*n)/3.*gamma*pow(3*alpha,n)*pow(H,2*(n-1));

	return X*(-(1-2*n)*(3-2*n)/6*gamma*pow(3*alpha,n)*pow(H,2*n)-H*H+3*Theta_Lambda)/T1;
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
	double H;

	H=Xdot/X;

	return H*H+(1.-2.*n)/6.*gamma*pow(3*alpha,n)*pow(H,2*n)-Theta_Lambda-Theta_m/(X*X*X);
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

double Theta_msolve(double X, double Y, double Xdot, double Ydot, double Theta_Lambda, double alpha, double gamma)
{
	double H;

	H=Xdot/X;

	return X*X*X*(H*H+(1-2*n)/6.*gamma*pow(3*alpha,n)*pow(H,2*n)-Theta_Lambda);
}





double dsolve(double Theta_m, double Theta_Lambda, double alpha, double gamma, double X_max)
{
	double X, Xdot, X1, Xdot1, Y, Ydot, Y1, Ydot1, zeta, mu_th, dL;
	double tau, tau_check, x, y, z, deviation;
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
		if(X>X_infty || X<X_max || i>max_iter || isnan(X))
		{
			printf("pop \n");
			break;
		}

		deviation=epsilon(X,Y,Xdot,Ydot, Theta_m, Theta_Lambda, alpha, gamma);
		if(fabs(deviation)>EPS)
		{
			printf("fuck \n");
			break;
		}
		/*
		if(i>max_iter)
		{
			break;
		}
		if(isnan(X))
		{
			break;
		}
		*/

                
                tau_check=(tau/h)/filter;
                if(floor(tau_check)==tau_check)
                {
			//deviation=epsilon(X,Y,Xdot,Ydot, Theta_m, Theta_Lambda, alpha, gamma);
			dL=zeta/X;
			mu_th=52.38560626+ 5*log10(dL); //Hub=25+5*log10(speed of light in km/s)
			//printf("%.10f %.10f %.10f %.10f %.10f %d \n",X,tau*H0, deviation, mu_th, 1./X-1.,i);
                }
		
	}	
	
	dL=zeta/X;
        mu_th=52.38560626 + 5*log10(dL); //Hub=52.38560626=25+5*log10(speed of light in km/s)
	return mu_th;
}

int main()
{
	double Theta_msub, Theta_Lambdasub, alphasub, gammasub, dL, mu, z, q;

	double T1, T2;
	
	//Parameter setup
	Theta_Lambdasub=Theta_Lambda_top; //0.5*H0*H0;
	alphasub=1.0;
	gammasub=gamma_top; //-10*0.65*H0*H0;
	Theta_msub=Theta_msolve(1.0, 1.0, H0, 1.0, Theta_Lambdasub, alphasub, gammasub);
	
	


	/* //diagnostics for deceleration parameter
	T1=(1-2*n)*(3-2*n)*gammasub*pow(3*alphasub,n)*pow(H0,2*n)/6;
	T2=2+n*(1-2*n)/3*gammasub*pow(3*alphasub,n)*pow(H0,2*(n-1));
	q=-Xddot(X_init, Y_init, Xdot_init, Ydot_init,Theta_msub,Theta_Lambdasub,alphasub,gammasub)/H0/H0;
	printf("%f %f %f %f \n" ,q, T1, T2, Theta_msub);
	*/  

	//z=200.0; //Given a z, calculate mu
	FILE *fp;
	int i=0,j=0;
	double arr_z[sizeofz];
	double arr_mu_th[sizeofz];
	fp=fopen("C:\\Users\\anwen\\Downloads\\SCPUnion2_1_z.txt","r");

	while(fp!=NULL)
    {
        fscanf(fp,"%lf",&arr_z[i]);
        z=arr_z[i];
        mu=dsolve(Theta_msub, Theta_Lambdasub, alphasub, gammasub, 1.0/(1+z));
        arr_mu_th[j]=mu;

        printf("%lf      ",arr_z[i]);
        printf("%lf\n",arr_mu_th[j]);

        if((arr_z[i]<=0))
        {
            printf("error number\n");
            break;
        }
        i++;
        j++;
    }
    fclose(fp);
	
	//printf("%.10f %.10f \n",z,mu);

	return 0;
}
