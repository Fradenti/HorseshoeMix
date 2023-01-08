#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////
// provides the matrix of cluster indexes, parameters over the cols, simulations over the rows
// [[Rcpp::export]]
arma::mat PSM(arma::mat inds){
  int nsim = inds.n_rows;
  int n    = inds.n_cols;
  arma::mat PSM(n,n), D(n,n);
  D.eye(); PSM.zeros();

  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      arma::colvec Z = (inds.col(i)-inds.col(j));
      arma::uvec success = find(Z==0);
      PSM(i,j) = success.n_elem;
      PSM(j,i) = PSM(i,j);
    }
  }
  return(PSM/nsim + D);
}
//-----------------------------------------------------------------------
// Various functions for sampling cluster indicators, for different scenarios
// 1 - classic case
// [[Rcpp::export]]
arma::colvec update_zeta_cpp(arma::colvec beta,
                         arma::colvec lambda2,
                         arma::colvec pi,
                         double tau2,
                         double sigma2,
                         int L, int p,
                         arma::colvec posslab){

  arma::colvec newz(p);
  arma::colvec sds = sqrt(lambda2 * tau2 * sigma2);
  //Rcout << sds << "\n";
  arma::colvec means(L,arma::fill::zeros);
  //Rcout << means << "\n";
  arma::colvec temp_prob(L);

  for(int pp=0;pp<p;pp++){

    arma::colvec Beta(L);
    Beta.fill(beta[pp]);

    temp_prob = arma::log_normpdf(Beta, means, sds) ;
    temp_prob = log(pi) + temp_prob;
    temp_prob = exp( temp_prob - max(temp_prob));
    newz[pp]  = RcppArmadillo::sample(posslab, 1, TRUE, temp_prob)[0];
    }

  return(newz);
}

// 2 - no sigma in beta specification
// [[Rcpp::export]]
arma::colvec update_zeta_cpp_nosigma(arma::colvec beta,
                             arma::colvec lambda2,
                             arma::colvec pi,
                             double tau2,
                             int L, int p,
                             arma::colvec posslab){

  arma::colvec newz(p);
  arma::colvec sds = sqrt(lambda2 * tau2);
  //Rcout << sds << "\n";
  arma::colvec means(L,arma::fill::zeros);
  //Rcout << means << "\n";
  arma::colvec temp_prob(L);

  for(int pp=0;pp<p;pp++){

    arma::colvec Beta(L);
    Beta.fill(beta[pp]);

    temp_prob = arma::log_normpdf(Beta, means, sds) ;
    temp_prob = log(pi) + temp_prob;
    temp_prob = exp( temp_prob - max(temp_prob));
    newz[pp]  = RcppArmadillo::sample(posslab, 1, TRUE, temp_prob)[0];
  }

  return(newz);
}

// 3 -  sampler with beta marginalized out
// [[Rcpp::export]]
arma::colvec update_zeta_cpp_marginalized_MEAN(arma::colvec Y,
                             arma::colvec lambda2,
                             arma::colvec pi,
                             double tau2, double sigma2,
                             int L, int n, arma::colvec posslab){

  arma::colvec newz(n);
  arma::colvec sds = sqrt( ( lambda2 * tau2 + 1 ) * sigma2 );
  arma::colvec means(L,arma::fill::zeros);

  for(int pp=0;pp<n;pp++){
  arma::colvec temp_prob(L);
  arma::colvec YY(L);

    YY.fill(Y[pp]);

    temp_prob =   arma::log_normpdf(YY, means, sds ) ;
    //Rcout<< "\n"<< temp_prob<<"\n";
    temp_prob = log(pi) + temp_prob;
    temp_prob = exp( temp_prob - max(temp_prob));
    newz[pp]  = RcppArmadillo::sample(posslab, 1, TRUE, temp_prob)[0];
  }

  return(newz);
}

//-----------------------------------------------------------------------
// Function for sampling Lambdas

// [[Rcpp::export]]
arma::colvec Slice_lambda_cpp(
                          arma::colvec beta,
                          arma::colvec prevlambda2,
                          arma::colvec zeta,
                          arma::colvec nl,
                          double sigma2,
                          double tau2,
                          int L){

  arma::colvec newlambda2(L);

  for( int l=0; l<L; l++ ){

    if(nl[l]==0){
      arma::colvec X = arma::randn<arma::colvec>(2);
      newlambda2[l]  = pow((X[1] / X[0]),2);
    }else{

      double t       = 1/prevlambda2[l];
      double u       = arma::randu(1)[1] * (1/(1+t));
      double LOWu    = 1/u -1;
      arma::uvec ind = find(zeta==l+1);

      double rate       = accu(pow( beta.elem(ind), 2)) / (2*sigma2*tau2);
      double F_uppBound     = R::pgamma(LOWu,(nl[l] + 1)/2, 1/rate,
                                        1,0);

    if(F_uppBound<.0001){F_uppBound     = 0.0001; }
      double p_g     = arma::randu(1)[0] *F_uppBound;

      double Quant      = R::qgamma(p_g,(nl[l] + 1)/2., 1./rate, 1, 0);
      newlambda2[l]  = 1/Quant;

    }
    }
    return(newlambda2);
  }


// [[Rcpp::export]]
arma::colvec Slice_lambda_scaled_cpp(
    arma::colvec beta,
    arma::colvec prevlambda2,
    arma::colvec zeta,
    arma::colvec nl,
    double sigma2,
    double tau2,
    double gamma2_lambda,
    int L){

  arma::colvec newlambda2(L);

  for( int l=0; l<L; l++ ){

    if(nl[l]==0){
      arma::colvec X = arma::randn<arma::colvec>(2);
      newlambda2[l]  = gamma2_lambda * pow((X[1] / X[0]),2);
    }else{

      double t       = gamma2_lambda/prevlambda2[l];
      double u       = arma::randu(1)[1] * (1/(1+t));
      double LOWu    = 1/u -1;
      arma::uvec ind = find(zeta==l+1);

      double rate       = accu(pow( beta.elem(ind), 2)) / (2*sigma2*tau2*gamma2_lambda);
      double F_uppBound     = R::pgamma(LOWu,(nl[l] + 1)/2, 1/rate,
                                        1,0);

      if(F_uppBound<.0001){F_uppBound     = 0.0001; }
      double p_g     = arma::randu(1)[0] *F_uppBound;

      double Quant      = R::qgamma(p_g,(nl[l] + 1)/2., 1./rate, 1, 0);
      newlambda2[l]  = gamma2_lambda/Quant;

    }
  }
  return(newlambda2);
}



//-----------------------------------------------------------------------

// [[Rcpp::export]]
double rivngamma_cpp(double a, double b){
  return 1/R::rgamma(a,1/b);
}

//-----------------------------------------------------------------------
// BNP

// [[Rcpp::export]]
arma::colvec UPD_Sticks_DP_cpp(arma::colvec Z,
                               int L,
                               double alphaDP){

  arma::colvec Sticks(L);

  for(int l=0; l<(L); l++){ //   for(int l=0; l<(L_new-1); l++){
    Sticks[l] = R::rbeta(1 + accu( Z == (l+1)),
                         alphaDP + accu(Z > (l+1)) );
  }

  return Sticks;
}


// [[Rcpp::export]]
arma::colvec UPD_Sticks_PY_cpp(arma::colvec Z,
                               int L,
                               double dPY, double alphaPY){

  arma::colvec Sticks(L);

  for(int l=0; l<(L); l++){ //   for(int l=0; l<(L_new-1); l++){
    Sticks[l] = R::rbeta(1 - dPY + accu( Z == (l+1)),
                         alphaPY + (l+1)*dPY + accu(Z > (l+1)) );
  }

  return Sticks;
}


// [[Rcpp::export]]
arma::colvec StickBreaker_cpp(arma::colvec V) {
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2   = exp(arma::cumsum(log(mV)));
  mV2   = arma::shift(mV2,+1);
  mV2(0) = 1;
  return(V%mV2);
}

