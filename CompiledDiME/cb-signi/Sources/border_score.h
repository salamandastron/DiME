
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_erf.h>
#include "./standard_package/standard_include.cpp"

/*
double f(double x) {
	
	double a=sqrt(6.2830);
	return 1/a*exp(-(x*x)/2);

}

*/




double order_cum(int n, int r, const double & x) {
	
	if(x>1)
		return 1;
	return gsl_cdf_binomial_Q(r-1, x, n);

}


/*


double order_pdf(int n, int r, const double & x) {
	
	double step=1e-7;
	return (order_cum(n, r, x+step) - order_cum(n, r, x))/step;

}
*/



double integrand(int n, int b, double x, double z) {

		
	if(x==0)
		return 0;

	
	double m= (b-1) * x/2.;
	double s= sqrt((b-1) /12.)*x;
	double t= ((z-x)-m)/s;
	
	
	return order_cum(n, b, x) / sqrt(2. * M_PI) * exp(-t*t/2.) *(z/(s*x));


}


double ordintegral(double a1, double a2, int steps, int n, int b, double z) {
	
	
	double step=(a2-a1)/steps;
	
	double integ=0;
	
	for(int i=0; i<=steps; i++) {
		
		if(i==0 || i==steps)
			integ+= 0.5 *integrand(n, b, i*step+a1, z)*step;
		else
			integ+= integrand(n, b, i*step+a1, z)*step;
	
	}
	
	
	return integ;	


}


double errf_from_minus_inf_to_x(double x, int b, double zeta) {
	
	double m= (b-1) * x/2.;
	double s= sqrt((b-1) /12.)*x;
	double t= (zeta-x-m)/s;
	
	return gsl_sf_erf_Q(-t);


}





double order_cum_bisection (double (*f)(int, int, const double &), double zerof, double a, double b, double precision, double & lo, double & hi, int n, int r) {
	
	
	
	if((f(n, r, a)-zerof)*(f(n, r, b)-zerof)>0) {
		cerr<<"wrong interval"<<endl;
		return a;
	}
	
	if ((f(n, r, a)-zerof) <= 0) {
		lo = a;
		hi = b;
	}
	else {
		lo = b;
		hi = a;
	}
	
	double mid = lo + (hi-lo)/2;
	
	while ((mid != lo) && (mid != hi) ) {
		
		
		double fmid=f(n, r, mid);
		
		if ((fmid-zerof) <= 0)
			lo = mid;
		else
			hi = mid;
		
		mid = lo + (hi-lo)/2;
		
		if(fabs(fmid)<precision) if(fabs(f(n, r, lo)-zerof)<precision && fabs(f(n, r, hi)-zerof)<precision)
			break;
		
		
	}
	

	return mid;
	
	
	
	
	
}


int integral_limits(int n, int b, double & lsx, double & ldx) {
	
	double lo, hi;
	order_cum_bisection(order_cum, 1-1e-7, 0, 1, 1e-8, lo, hi, n, b);
	ldx=lo;
	order_cum_bisection(order_cum, 1e-7, 0, 1, 1e-8, lo, hi, n, b);
	lsx=hi;

	return 0;

}



double cumformula_gaussian(double z_z, int n, double b, double r0) {
	
	
	
	double zeta= (z_z - b * r0) / (1 - r0);
	if(zeta<0)
		return 0;
	
	if(zeta>b) {
	
		cerr<<"ERROR in the border score!"<<endl;
		return 1;
		
	}
	
	if(b==1)
		return order_cum(n, b, zeta);
	if(b==2) {
	
		if(zeta<1)
			return 1 - 2 *pow((1-zeta/2.), n) + pow((1-zeta), n);
		else
			return 1 - 2 *pow((1-zeta/2.), n);
			
	
	}
		
	double lsx, ldx;
	integral_limits(n, b, lsx, ldx);

	double Icontr=ordintegral(lsx, ldx, 100, n, b, zeta);

	return  min(1., Icontr+ errf_from_minus_inf_to_x(ldx, b, zeta));
	
	
}



