//Because I don't create a package including these functions, I have to combine the fnb function in this file

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm> 
#include <cmath> 
#define BIG 1e20
#define MAXIT 50

using namespace Rcpp;

List pfn(NumericMatrix x, NumericVector y, double tau, double Mm_factor, int max_bad_fixup, double eps);
NumericMatrix sub_matrix(NumericMatrix mat, IntegerVector samp, bool row_or_col);
NumericVector sub_vector(NumericVector vec, IntegerVector samp);
NumericMatrix cross_chol_inv(NumericMatrix x);
NumericMatrix times_sqrt(NumericMatrix x, NumericMatrix xxinv, NumericVector p);
IntegerVector trueVector(LogicalVector x);
NumericMatrix rbind(NumericMatrix mat, NumericVector vec);
NumericVector cbind(NumericVector vec, double val);
NumericVector glob_transform(NumericMatrix mat, LogicalVector vec);
LogicalVector combine_logic(LogicalVector a, LogicalVector b);
int sum_or(LogicalVector a, LogicalVector b);
NumericVector mat_times_vec(NumericMatrix x, NumericVector p);

//fnb
void stepy(arma::vec &b, arma::vec &d, arma::mat &a, arma::mat &ada);
arma::vec lpfnb(int p, int n, arma::mat &a,arma::vec &b, arma::vec &c, arma::vec &x, double beta, double eps);
void rqfnb(NumericMatrix x, NumericVector y, double tau, double beta, double eps);

//Call R function
NumericVector quantile(NumericVector x, NumericVector probs);



//[[Rcpp::export]]
List pfn(NumericMatrix x, NumericVector y, double tau = 0.5,  double Mm_factor = 0.8, int max_bad_fixup = 3, double eps  = 1e-06){

	int n = y.length();
	if (x.nrow() != n) 
		stop("X and y don't match n");
	if (tau<0 || tau>1)
		stop("tau outside (0,1)");
	int p = x.ncol();
	int m = round( pow( (p+1.0)*n, 2.0/3 ) ); //Rcout<< "m= " <<m<<std::endl;
	bool not_optimal = true;
	NumericVector b(p);

	while (not_optimal){

		if ( m>=n ){
			return rqfnb(x,y,tau,0.99995,eps);//call fnb
		}
		else{
				IntegerVector tmp = seq_len(n); 
				IntegerVector s = RcppArmadillo::sample(tmp,m,false, NumericVector::create()); 
				NumericMatrix xx = sub_matrix(x,s-1,0);// the index in C++ starts at 0
				NumericVector yy = sub_vector(y,s-1);
				List z = rqfnb(xx,yy,tau,0.99995,eps);
				NumericMatrix xxinv = cross_chol_inv(xx);
				NumericMatrix band = times_sqrt(x,xxinv,NumericVector(p,1.0));
				NumericVector r = y-mat_times_vec(x,z["coefficients"]);
				double M = Mm_factor*m;
				double lo_q = std::max(1.0/n,tau-M/(2.0*n));
				double hi_q = std::min(tau+M/(2.0*n),(n-1.0)/n);
				NumericVector kappa = quantile(r/pmax(eps,band),NumericVector::create(lo_q,hi_q));
				LogicalVector su = (r > ( band*kappa(1) ) );
				LogicalVector sl = (r < ( band*kappa(0) ) );
				int bad_fixup = 0;
				while (not_optimal && ( bad_fixup < max_bad_fixup ) ){
					xx = sub_matrix(x,trueVector( (!su) & (!sl) ), 0);
					yy = sub_vector(y,trueVector( (!su) & (!sl) ) );
					if ( is_true(any(su)) ){
						NumericVector ghib_x = glob_transform(x,su);
						double ghib_y = sum( sub_vector(y,trueVector(su)) );
						xx = rbind(xx,ghib_x);
						yy = cbind(yy,ghib_y);
					}
					if ( is_true(any(sl)) ){
						NumericVector glob_x = glob_transform(x,sl);
						double glob_y = sum( sub_vector(y,trueVector(sl)) );
						xx = rbind(xx,glob_x);
						yy = cbind(yy,glob_y);
					}
					z = rqfnb(xx,yy,tau,0.99995,eps);
					b = z["coefficients"];
					r = y-mat_times_vec(x,z["coefficients"]);
					LogicalVector su_bad = (( r<0 ) & su);
					LogicalVector sl_bad = (( r>0 ) & sl);
					if ( is_true(any(combine_logic(su_bad,sl_bad))) ) {
						if ( sum_or(su_bad,sl_bad) > 0.1*M){
							warning("Too many fixups: doubling m");
							m = 2*m;
							break;
						}	
						su = su & (!su_bad);
						sl = sl & (!sl_bad);
						bad_fixup++;
					}
					else
						not_optimal = false;
				}
		}
	}
	NumericVector coefficients(b);
	coefficients.attr("names") = colnames(x);
	return List::create(Named("coefficients") = coefficients, Named("tau") = tau, Named("residuals") = y-mat_times_vec(x,coefficients));

}


NumericMatrix cross_chol_inv(NumericMatrix x){

	arma::mat x_arma(x.begin(),x.nrow(),x.ncol(),false);

	arma::mat tmp = arma::trans(x_arma)*x_arma;

	arma::mat rt = arma::inv( arma::chol( tmp ) );

	return as<NumericMatrix>(wrap(rt));
}

NumericMatrix times_sqrt(NumericMatrix x, NumericMatrix xxinv, NumericVector p){

	arma::mat x_arma(x.begin(),x.nrow(),x.ncol(),false);
	arma::mat xxinv_arma(xxinv.begin(),xxinv.nrow(),xxinv.ncol(),false);
	arma::mat p_arma(p.begin(),p.length(),1,false);

	arma::mat tmp = arma::pow((x_arma)*xxinv_arma,2);
	arma::mat rt = arma::sqrt(tmp*p_arma);

	return as<NumericMatrix>(wrap(rt));
}

//copy the sub-matrix into new matrix
NumericMatrix sub_matrix(NumericMatrix mat, IntegerVector samp, bool row_or_col){
	
	if (!row_or_col){

		int row = samp.length();
		int col = mat.ncol();
		IntegerVector sub_vec( col );
		NumericMatrix rt( row,col );
	
		// need modification by Range()
		for (int i = 0; i< row; i++){
			rt(i,_) = mat(samp[i],_);
		}
		return rt;
	}
	else{

		int row = mat.nrow();
		int col = samp.length();
		IntegerVector sub_vec( row );
		NumericMatrix rt( row,col );
	
		for (int i = 0; i< col; i++){
			rt(_,i) = mat(_,samp[i]);
		}
		return rt;
	}

}

//copy the sub-vector into new vector
NumericVector sub_vector(NumericVector vec, IntegerVector samp){
	
	int len = samp.length();
	NumericVector rt(len);
	for (int i=0; i<len; i++){
		rt(i) = vec(samp[i]);
	}
	return rt;
}

NumericVector quantile(NumericVector x, NumericVector probs) {

	Environment stats("package:stats");
	Function quantile = stats["quantile"];
	int npr = probs.size();
	NumericVector ans(npr);
	for(int i=0; i<npr; i++){
		ans[i] = as<double>(quantile(x, probs[i]));
  	}
  	return ans;
}

IntegerVector trueVector(LogicalVector x) {

	std::vector<int> count;
	for (int i =0; i<x.length(); i++){
		if (x[i] == true)
			count.push_back(i);
	}

	return as<IntegerVector>(wrap(count));
}

NumericMatrix rbind(NumericMatrix mat, NumericVector vec){

	if (mat.ncol() != vec.length())
		stop("matrix and vector don't match col");
	int new_nrow = mat.nrow()+1;
	NumericMatrix rt(new_nrow,mat.ncol());
	
	for (int i = 0; i< new_nrow-1; i++){
		rt(i,_) = mat(i,_);
	}	
	rt(new_nrow-1,_) = vec;
	return rt;
}

NumericVector cbind(NumericVector vec, double val){

	int new_len = vec.length()+1;
	NumericVector rt(new_len);
	for (int i = 0; i< new_len-1; i++){
		rt[i] = vec[i];
	}	
	rt[new_len-1] = val;
	return rt;
}

NumericVector glob_transform(NumericMatrix mat, LogicalVector vec){

	NumericMatrix sub = sub_matrix(mat,trueVector(vec),0);
	NumericVector rt(sub.ncol());

	for (int i =0; i<rt.length(); i++){
		for (int j = 0; j < sub.nrow(); j++)
			rt[i] += sub(j,i);
	}
	return rt;
}

LogicalVector combine_logic(LogicalVector a, LogicalVector b){

	int na = a.length(), nb = b.length();
	LogicalVector rt(na+nb);
	for (int i = 0; i<na ; i++)
		rt[i] = a[i];
	for (int j = 0; j<nb ; j++)
		rt[na+j] = b[j];
	return rt;
}

int sum_or(LogicalVector a, LogicalVector b){

	LogicalVector or_gate = a|b;
	int count = 0;
	for (int i =0; i<or_gate.length(); i++)
		if (or_gate[i] == true)
			count++;
	return count;
}

NumericVector mat_times_vec(NumericMatrix x, NumericVector p){

	arma::mat x_arma(x.begin(),x.nrow(),x.ncol(),false);
	arma::mat p_arma(p.begin(),p.length(),1,false);
	return as<NumericVector>(wrap(x_arma*p_arma));
}


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
