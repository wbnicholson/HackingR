

#include <RcppArmadillo.h>
#include <math.h> 

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//Elementwise Soft-thresholding Function
// [[Rcpp::export]]
double ST1(double z,double gam){
double sparse=0;
if(z>0 && gam<fabs(z)) return(z-gam);

if(z<0 && gam<fabs(z)) return(z+gam);
if(gam>=fabs(z)) return(sparse);
else return(0);
}

//Soft-thresholding by column
// [[Rcpp::export]]
colvec ST3(colvec z ,double gam)
{
int n=z.size();
colvec z1(n);
for( int i=0; i<n;++i)
{
double z11=z(i);
z1(i)=ST1(z11,gam);
}
return(z1);
}

//Subsetting to calculate partial residual. 
//Use of the Rcpp Integer vector causes problems in parallelization
uvec ind(int n2,int m){
IntegerVector subs(n2);
for(int i =0 ; i<n2;++i)
{
subs(i)=i;
 }
subs.erase(m);
return(as<uvec>(subs));
}

//The actual lasso function 
// [[Rcpp::export]]
mat lassocore(mat Y,mat Z, mat B,mat BOLD,mat R, const double gam,const double lassothresh,colvec ZN,colvec Z2, colvec G1, NumericVector G, uvec m, int k2,int n2){
double thresh=10*lassothresh;

while(thresh>lassothresh)
{

  for(int i = 0; i < k2; ++i)
	{
	m=ind(k2,i);	       
        		
	R=Y-B.cols(m)*Z.rows(m);
	Z2=trans(Z.row(i));
     	G1=R*Z2;
	G=as<NumericVector>(wrap(G1));
	B.col(i)=ST3(G,gam)/ZN[i];

	  }	

 mat thresh1=abs((B-BOLD)/(ones(n2,k2)+abs(BOLD)));
 thresh=norm(thresh1,"inf");
 BOLD=B;
 }

return(B);
}
//Loops through grid penalty values
// [[Rcpp::export]]
List gamloop(List beta, NumericMatrix y,NumericMatrix z, NumericVector gamm, const double lassothresh,NumericVector YMean, NumericVector ZMean,NumericVector znorm2,NumericMatrix BFoo){
colvec YMean2=YMean;
colvec ZMean2=ZMean;
int i;
NumericMatrix BOLD=BFoo;

double gammgrid[gamm.size()];
int gran2=gamm.size();

int n=y.nrow(),k=y.ncol();
arma::mat Y(y.begin(),n,k,false);
int n1=z.nrow(),k1=z.ncol();
arma::mat Z(z.begin(),n1,k1,false);
int n2=BFoo.nrow(),k2=BFoo.ncol();
arma::mat B1(BFoo.begin(),n2,k2,false);
 mat b2=B1;
 mat B1F2;
 mat R=zeros<mat>(n,k);
colvec ZN = as<colvec>(znorm2);
 colvec nu=zeros<colvec>(n2);
 colvec Z2=trans(Z.row(1));
colvec G1=R*Z2;
 NumericVector G=as<NumericVector>(wrap(G1));
 uvec m=ind(k2,1);

 double gam =0;
for(int j=0;j<gran2;++j)
{
gammgrid[j]=gamm[j];
}

   for (i=0; i<gamm.size();++i) {
	gam=gammgrid[i];
	mat B1F2(as<NumericMatrix>(beta[i]).begin(),n2,k2,false);
 	B1 = lassocore(Y,Z,B1F2,b2,R,gam, lassothresh,ZN,Z2,G1,G,m,k2,n2); 
 	nu = YMean2 - B1 *ZMean2;
	beta[i] = mat(join_horiz(nu, B1)); 
  
}
    return(wrap(beta));
}


