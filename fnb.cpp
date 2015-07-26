// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm> 
#include <cmath> 
#define BIG 1e20
#define MAXIT 50

using namespace Rcpp;


void stepy(arma::vec &b, arma::vec &d, arma::mat &a, arma::mat &ada){

	arma::mat tmp = arma::trans(a);
	for (int i =0; i<(int)tmp.n_rows; i++){
		tmp.row(i) *= d(i);
	}

	ada = a*tmp;//dsyr('U',p,d(i),a(1,i),1,ada,p)

	b = arma::solve(ada,b);//dposv('U',p,1,ada,p,b,p,info)

	ada = arma::chol(ada);
	/* what dposv does; possibly faster
	ada = chol(ada);
	b = arma::solve(trans(trimatu(ada))*trimatu(ada),b);//dposv('U',p,1,ada,p,b,p,info)
	*/
}

arma::vec lpfnb(int p, int n, arma::mat &a,arma::vec &b, arma::vec &c, arma::vec &x, double beta, double eps){

	arma::vec nit(3,arma::fill::zeros);
	nit(2) = n;
	arma::vec d(n,arma::fill::ones);//rep(1,n)
	arma::vec y = a*c;//dgemv('N',p,n,one,a,p,c,1,zero,y,1)
					  //actually it's gonna be coefficients
	arma::mat ada(p,p,arma::fill::zeros);

	stepy(y,d,a,ada);

	arma::vec s = -arma::trans(a)*y+c;//dcopy(n,c,1,s,1) + dgemv('T',p,n,mone,a,p,y,1,one,s,1)

	arma::vec u(n,arma::fill::ones);
	arma::vec z(n,arma::fill::zeros);
	arma::vec w(n,arma::fill::zeros);
	for (int i = 0; i < n; i++){	
		if (std::abs(s(i)) < eps){
			z(i) = std::max(s(i),0.0)+eps;
			w(i) = std::max(-s(i),0.0)+eps;
		}
		else{
			z(i) = std::max(s(i),0.0);
			w(i) = std::max(-s(i),0.0);
		}
		s(i) = 1 - x(i);
	}

	double gap = arma::dot(z,x)+arma::dot(w,s);

	arma::vec dx(n,arma::fill::zeros);
	arma::vec ds(n,arma::fill::zeros);
	arma::vec dz(n,arma::fill::zeros);
	arma::vec dy(n,arma::fill::zeros);
	arma::vec rhs(n,arma::fill::zeros);
	arma::vec dw(n,arma::fill::zeros);
	arma::vec dr(n,arma::fill::zeros);

	while ((gap > eps) && (nit(0) < MAXIT) ){

		nit(0) ++;
		for (int i = 0; i < n; i++){
			d(i) = 1.0/(z(i)/x(i)+w(i)/s(i));
			ds(i) = z(i)-w(i);
			dz(i) = d(i)*ds(i);
		}

		dy = -a*x+b;//dcopy(p,b,1,dy,1) + dgemv('N',p,n,mone,a,p,x,1,one,dy,1)
		dy = a*dz+dy;//dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
		rhs = dy;//dcopy(p,dy,1,rhs,1) 

		stepy(dy,d,a,ada);
	

		ds = arma::trans(a)*dy-ds;//call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)

		double deltap = BIG;
		double deltad = BIG;
		for (int i = 0; i<n; i++){
			dx(i)=d(i)*ds(i);
			ds(i)=-dx(i);
			dz(i)=-z(i)*(dx(i)/x(i) + 1.0);
			dw(i)=-w(i)*(ds(i)/s(i) + 1.0);
			if(dx(i) < 0)
	      		deltap=std::min(deltap,-x(i)/dx(i));
	      	if(ds(i) < 0)
	      		deltap=std::min(deltap,-s(i)/ds(i));
	      	if(dz(i) < 0)
	      		deltad=std::min(deltad,-z(i)/dz(i));
	      	if(dw(i) < 0)
	      		deltad=std::min(deltad,-w(i)/dw(i));
		}
		deltap = std::min(beta*deltap,1.0);
		deltad = std::min(beta*deltad,1.0);

		if (std::min(deltap,deltad) < 1){
			nit(1) ++;
			double mu = arma::dot(z,x)+arma::dot(w,s);
			double g = mu + deltap*arma::dot(dx,z) + deltad*arma::dot(dz,x) + deltap*deltad*arma::dot(dx,dz)
			+ deltap*arma::dot(ds,w) + deltad*arma::dot(dw,s) + deltap*deltad*arma::dot(ds,dw);//gap is mu
			mu = mu*std::pow(g/mu,3)/(2.0*n);
			for (int i =0; i < n; i++){
				dr(i) = d(i)*(mu*(1/s(i)-1/x(i))+ dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i));
			}
			std::swap(rhs,dy); 
			dy = a*dr + dy; //dgemv('N',p,n,one,a,p,dr,1,one,dy,1)
			dy = arma::solve(arma::trans(ada)*ada,dy); //dpotrs('U',p,1,ada,p,dy,p,info)
			u = arma::trans(a)*dy; //dgemv('T',p,n,one,a,p,dy,1,zero,u,1)
	      	deltap=BIG;
	      	deltad=BIG;
	      	for (int i = 0; i<n; i++){
	      		double dxdz = dx(i)*dz(i);
	      		double dsdw = ds(i)*dw(i);
	      		dx(i)= d(i)*(u(i)-z(i)+w(i))-dr(i);
	     		ds(i)= -dx(i);
	      		dz(i)= -z(i)+(mu - z(i)*dx(i) - dxdz)/x(i);
	      		dw(i)= -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i);
	      		if(dx(i) < 0)
	      			deltap=std::min(deltap,-x(i)/dx(i));
	      		if(ds(i) < 0)
	      			deltap=std::min(deltap,-s(i)/ds(i));
	      		if(dz(i) < 0)
	      			deltad=std::min(deltad,-z(i)/dz(i));
	      		if(dw(i) < 0)
	      			deltad=std::min(deltad,-w(i)/dw(i));
	      	}

	      	deltap=std::min(beta*deltap,1.0);
	      	deltad=std::min(beta*deltad,1.0);
	    }
	    x = deltap*dx + x; //daxpy(n,deltap,dx,1,x,1)
	    s = deltap*ds + s; //daxpy(n,deltap,ds,1,s,1)
	    y = deltad*dy + y; //daxpy(p,deltad,dy,1,y,1)
	    z = deltad*dz + z; //daxpy(n,deltad,dz,1,z,1)
		w = deltad*dw + w; //daxpy(n,deltad,dw,1,w,1)
	    gap = arma::dot(z,x)+arma::dot(w,s);
	}
/*
	z = -w+z; //daxpy(n,mone,w,1,z,1)
	std::swap(x,z);
*/
	return y;
}

//[[Rcpp::export]]
void rqfnb(NumericMatrix x, NumericVector y, double tau, double beta, double eps){

	int n = y.length();
	int p = x.ncol();
	if (x.nrow() != n) 
		stop("X and y don't match n");
	if (tau<eps || tau>(1-eps))
		stop("No parametric Frisch-Newton method. Set tau in (0,1)");

	arma::mat x_arma(x.begin(),n,p,false);
	arma::vec y_arma(y.begin(),n,false);
	arma::vec rhs(p);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < p; j++){
			rhs(j)+=x(i,j);
		}
	}
	rhs = rhs*(1-tau);
	arma::vec wn(n);
	wn.fill(1-tau);
	arma::mat a = arma::trans(x_arma);
	arma::vec c = -y_arma;
	arma::vec coef = lpfnb(p, n ,a, rhs, c, wn, beta, eps);

	NumericVector coefficients = as<NumericVector>(wrap(-coef));
	NumericVector res = as<NumericVector>(wrap(y_arma+x_arma*coef));

	coefficients.attr("names") = colnames(x);

	return List::create(Named("coefficients") = coefficients, Named("tau") = tau, Named("residuals") = res);
}
