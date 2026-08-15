
#define TMB_LIB_INIT R_init_phylosem
#include <TMB.hpp>

// SparseMatrix for Ornstein-Uhlenbeck network correlations
//template<class Type>
//Eigen::SparseMatrix<Type> Q_network( Type log_alpha,
//                                     int n_s,
//                                     vector<int> parent_s,
//                                     vector<int> child_s,
//                                     vector<Type> dist_s ){
//
//  Eigen::SparseMatrix<Type> Q( n_s, n_s );
//  Type alpha = exp( log_alpha );
//  for(int s=0; s<n_s; s++){
//    Q.coeffRef( s, s ) = Type(1.0);
//  }
//  for(int s=1; s<parent_s.size(); s++){
//    if( exp(-dist_s(s))!=0 ){
//      Q.coeffRef( parent_s(s), child_s(s) ) = -exp(-alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), parent_s(s) ) = Q.coeffRef( parent_s(s), child_s(s) );
//      Q.coeffRef( parent_s(s), parent_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), child_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//    }
//  }
//  return Q;
//}

// ICAR for BM with sum-to-zero constraint
//  See: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-09 -- ICAR specification

// Precision of evolutionary covariance
//template<class Type>
//Eigen::SparseMatrix<Type> Q_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 int n_vars ){
//
//  // Define temporary objects
//  Eigen::SparseMatrix<Type> Q_vv( n_vars, n_vars );
//  // SEM
//  Eigen::SparseMatrix<Type> Linv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Rho_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gamma_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gammainv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  for(int zI=0; zI<RAM.rows(); zI++){
//    if(RAM(zI,0)==1) Rho_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1);
//    if(RAM(zI,0)==2) Gamma_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1); // Cholesky of covariance, so -Inf to Inf;
//  }
//  Gammainv_vv = atomic::matinv( Gamma_vv );
//  Linv_vv = Gammainv_vv * ( I_vv - Rho_vv );
//  Q_vv = Linv_vv.transpose() * Linv_vv;
//  return Q_vv;
//}

// Evolutionary covariance
//template<class Type>
//matrix<Type> V_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 vector<Type> RAMstart,
//                                 int n_vars ){
//
//  // Define temporary objects
//  matrix<Type> V_vv(n_vars, n_vars);
//  // SEM
//  matrix<Type> L_vv(n_vars, n_vars);
//  matrix<Type> Rho_vv(n_vars, n_vars);
//  matrix<Type> Gamma_vv(n_vars, n_vars);
//  matrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  Type tmp;
//  for(int r=0; r<RAM.rows(); r++){
//    // Extract estimated or fixed value
//    if(RAM(r,3)>=1){
//      tmp = beta_z(RAM(r,3)-1);
//    }else{
//      tmp = RAMstart(r);
//    }
//    // Assign to proper matrix
//    if(RAM(r,0)==1){
//      Rho_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }else{
//      Gamma_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }
//  }
//  L_vv = I_vv - Rho_vv;
//  L_vv = atomic::matinv( L_vv );
//  L_vv = L_vv * Gamma_vv;
//  V_vv = L_vv * L_vv.transpose();
//  return V_vv;
//}

template<class Type>
vector<Type> get_nonzero_elements( Eigen::SparseMatrix<Type> M ) {
    std::vector<Type> values;

    for (int k = 0; k < M.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Type>::InnerIterator it(M, k); it; ++it) {
            values.push_back(it.value());
        }
    }

    return values;  // Return vector of nonzero values if needed
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //using namespace Eigen;
  using namespace density;

  // Data
  DATA_INTEGER( n_tip );
  DATA_IMATRIX( edge_ez );
  DATA_VECTOR( length_e );
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );
  DATA_INTEGER( estimate_ou );
  DATA_INTEGER( estimate_lambda );
  //DATA_INTEGER( estimate_kappa );
  DATA_VECTOR( height_v );
  DATA_MATRIX( y_ij );
  DATA_UPDATE( y_ij ); // Experiment with cAIC
  DATA_IVECTOR( v_i );
  DATA_IVECTOR( familycode_j );
  DATA_IVECTOR( linkcode_j );
  DATA_IVECTOR( sigmastart_j );
  DATA_IVECTOR( group_j );
  DATA_STRING( experiments_type );

  // Parameters
  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_z );
  PARAMETER( lnalpha );
  PARAMETER( logitlambda );
  PARAMETER( lnkappa );
  PARAMETER_MATRIX( x_vj );
  PARAMETER_VECTOR( xbar_j );

  // Indices
  int n_e = edge_ez.rows();   // edges
  int n_v = x_vj.rows();      // vertices = edges + 1
  int n_j = x_vj.cols();      // variables
  int n_i = y_ij.rows();      // data
  int n_g = group_j.maxCoeff();
  // n_p = number of populations ... optional
  // n_k = number of experimental observations ... optional
  int vparent, vchild;
  int vroot = n_tip;          // vroot = n_tip+1 - 1, where latter -1 is converting from R to CPP indexing

  // JNLL
  Type jnll = 0;
  vector<Type> jnll_v( n_v );         // Likelihood by phylogenetic edge
  matrix<Type> jnll_ij( n_i, n_j );   // Likelihood by sample in y_ij (excluding base level for grouped categories)
  matrix<Type> jnll_ig( n_i, n_g );   // Likelihood for base level per group of factor levels
  jnll_ij.setZero();
  jnll_v.setZero();
  jnll_ig.setZero();

  // Global vars
  Type alpha = exp( lnalpha );
  Type lambda = invlogit( logitlambda );
  Type kappa = exp( lnkappa );
  matrix<Type> yhat_ij( n_i, n_j );
  matrix<Type> mu_vj( n_v, n_j );
  vector<Type> rho_v( n_v );
  vector<Type> var_v( n_v );
  vector<Type> sigma_z = exp( lnsigma_z );
  matrix<Type> eps_vj( n_v, n_j );
  matrix<Type> sumpred_vg( n_v, n_g );
  matrix<Type> sumobs_ig( n_i, n_g );
  eps_vj.setZero();
  sumpred_vg.setZero();
  sumobs_ig.setZero();

  //// Assemble Evolutionary covariance
  //matrix<Type> Vtmp( n_j, n_j );
  //matrix<Type> V_jj( n_j, n_j );
  //matrix<Type> L_jj(n_j, n_j);
  //matrix<Type> Rho_jj(n_j, n_j);
  //matrix<Type> Gamma_jj(n_j, n_j);
  //matrix<Type> I_jj( n_j, n_j );
  //Rho_jj.setZero();
  //Gamma_jj.setZero();
  //I_jj.setIdentity();
  //Type tmp;
  //for(int r=0; r<RAM.rows(); r++){
  //  // Extract estimated or fixed value
  //  if(RAM(r,3)>=1){
  //    tmp = beta_z(RAM(r,3)-1);
  //  }else{
  //    tmp = RAMstart(r);
  //  }
  //  // Assign to proper matrix
  //  if(RAM(r,0)==1){
  //    Rho_jj( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
  //  }else{
  //    Gamma_jj( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
  //  }
  //}
  //L_jj = I_jj - Rho_jj;
  //L_jj = atomic::matinv( L_jj );
  //L_jj = L_jj * Gamma_jj;
  //V_jj = L_jj * L_jj.transpose();

  //// Assemble Evolutionary covariance
  Eigen::SparseMatrix<Type> Qtmp( n_j, n_j );
  Eigen::SparseMatrix<Type> Rho_jj(n_j, n_j);
  Eigen::SparseMatrix<Type> Gamma_jj(n_j, n_j);
  Eigen::SparseMatrix<Type> I_jj( n_j, n_j );
  I_jj.setIdentity();
  Type tmp;
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    // Assign to proper matrix
    if(RAM(r,0)==1){
      Rho_jj.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    }
    if(RAM(r,0)==2){
      Gamma_jj.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp; // Cholesky of covariance, so -Inf to Inf;
    }
    if(RAM(r,0)==3){
      // done by edge later
    }
  }
  Eigen::SparseMatrix<Type> V_jj = Gamma_jj.transpose() * Gamma_jj;
  matrix<Type> Vinv_jj = tmbutils::invertSparseMatrix( V_jj );
  Eigen::SparseMatrix<Type> Vinv2_jj = asSparseMatrix( Vinv_jj );

  // Distribution of OU evolution -- Root
  // Correlation between i and parent(i) as distance -> INF
  // Load in moderators for root
  for(int r=0; r<RAM.rows(); r++){
    if(RAM(r,0)==3){
      Rho_jj.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = x_vj( vroot, RAM(r,4)-1 );
    }
  }
  Eigen::SparseMatrix<Type> IminusRho_jj = I_jj - Rho_jj;
  Eigen::SparseMatrix<Type> Q_jj = IminusRho_jj.transpose() * Vinv2_jj * IminusRho_jj;

  if( estimate_ou==1 ){
    rho_v(vroot) = 0;
    // SD of Ornstein-Uhlenbeck process as distance -> INF
    var_v(vroot) = Type(1.0) / (2*alpha);
    // conditional probability
    for(int j=0; j<n_j; j++){
      eps_vj(vroot,j) = x_vj(vroot,j) - xbar_j(j);
    }
    //Vtmp = V_jj * var_v(vroot);
    //jnll_v(vroot) = MVNORM(Vtmp)( eps_vj.row(vroot) );
    Qtmp = Q_jj / var_v(vroot);
    jnll_v(vroot) = GMRF(Qtmp)( eps_vj.row(vroot) );
    // Optionally fix the root at the mean
  }else{
    rho_v(vroot) = NAN;
    var_v(vroot) = NAN;
  }
  // Distribution of OU evolution -- Edges
  for(int e=0; e<n_e; e++){ // PARALLEL_REGION
    vchild = edge_ez(e,1);
    vparent = edge_ez(e,0);

    // Load in moderators
    for(int r=0; r<RAM.rows(); r++){
      if(RAM(r,0)==3){
        Rho_jj.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 0.5 * ( x_vj( vchild, RAM(r,4)-1 ) + x_vj( vparent, RAM(r,4)-1 ) );
      }
    }
    IminusRho_jj = I_jj - Rho_jj;
    Q_jj = IminusRho_jj.transpose() * Vinv2_jj * IminusRho_jj;

    if( estimate_ou==1 ){
      // Correlation between i and parent(i)
      rho_v(vchild) = exp( -alpha * pow(length_e(e),kappa) );
      // SD of O-U process
      var_v(vchild) = Type(1.0)/(2*alpha) * (Type(1.0)-exp( -2 * alpha * pow(length_e(e),kappa) ));
    }else{
      rho_v(vchild) = Type(1.0);
      var_v(vchild) = pow(length_e(e),kappa);
    }
    // conditional probability
    if( estimate_lambda==1 ){
      if( vchild < n_tip ){
        //Vtmp = V_jj * ( lambda*var_v(vchild) + (1-lambda)*height_v(vchild) );
        Qtmp = Q_jj / ( lambda*var_v(vchild) + (1-lambda)*height_v(vchild) );
      }else{
        //Vtmp = V_jj * ( lambda*var_v(vchild) );
        Qtmp = Q_jj / ( lambda*var_v(vchild) );
      }
    }else{
      //Vtmp = V_jj * var_v(vchild);
      Qtmp = Q_jj / var_v(vchild);
    }
    //xtmp_j = (x_vj.row(vchild).array()-xbar_j) - rho_v(vchild)*(x_vj.row(vparent).array()-xbar_j);
    for(int j=0; j<n_j; j++){
      eps_vj(vchild,j) = (x_vj(vchild,j)-xbar_j(j)) - rho_v(vchild)*(x_vj(vparent,j)-xbar_j(j));
    }
    //jnll_v(vchild) = MVNORM(Vtmp)( eps_vj.row(vchild) );
    jnll_v(vchild) = GMRF(Qtmp)( eps_vj.row(vchild) );
  }
  jnll += jnll_v.sum();

  // Distribution for data
  for(int v=0; v<n_v; v++){
    for(int j=0; j<n_j; j++){
      // Link function
      if( linkcode_j(j)==0 ){
        // identity link
        mu_vj(v,j) = x_vj(v,j);
      }
      if( linkcode_j(j)==1 ){
        // log link
        mu_vj(v,j) = exp(x_vj(v,j));
      }
      if( linkcode_j(j)==2 ){
        // logit link
        mu_vj(v,j) = invlogit(x_vj(v,j));
      }
      if( linkcode_j(j)==3 ){
        // cloglog link
        mu_vj(v,j) = Type(1.0) - exp( -1.0 * exp(x_vj(v,j)) );
      }
      if( linkcode_j(j)==4 ){
        // multivariate-logit link
        mu_vj(v,j) = exp(x_vj(v,j));
        sumpred_vg(v,group_j(j)-1) += exp(x_vj(v,j));
      }
    }
    // Divide by total across levels for each group for categorical traits
    for(int j=0; j<n_j; j++){
      if( linkcode_j(j)==4 ){
        // multivariate-logit link
        mu_vj(v,j) /= 1.0 + sumpred_vg(v,group_j(j)-1) ;
      }
    }
  }

  // Distribution for data
  for(int i=0; i<n_i; i++){
  for(int j=0; j<n_j; j++){
    yhat_ij(i,j) = mu_vj(v_i(i),j);
    // Likelihood
    //if( familycode_j(j)==0 ){
    // familycode = 0 :  don't include likelihood
    //}
    if( familycode_j(j)==1 ){
      // familycode = 1 :  normal
      if(R_FINITE(asDouble(y_ij(i,j)))){
        jnll_ij(i,j) -= dnorm( y_ij(i,j), mu_vj(v_i(i),j), sigma_z(sigmastart_j(j)), true );
      }
    }
    if( familycode_j(j)==2 ){
      // familycode = 2 :  binomial
      if(R_FINITE(asDouble(y_ij(i,j)))){
        jnll_ij(i,j) -= dbinom( y_ij(i,j), Type(1.0), mu_vj(v_i(i),j), true );
      }
    }
    if( familycode_j(j)==3 ){
      // familycode = 3 :  Poisson
      if(R_FINITE(asDouble(y_ij(i,j)))){
        jnll_ij(i,j) -= dpois( y_ij(i,j), mu_vj(v_i(i),j), true );
      }
    }
    if( familycode_j(j)==4 ){
      // familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
      if(R_FINITE(asDouble(y_ij(i,j)))){
        jnll_ij(i,j) -= dgamma( y_ij(i,j), pow(sigma_z(sigmastart_j(j)),-2), mu_vj(v_i(i),j)*pow(sigma_z(sigmastart_j(j)),2), true );
      }
    }
    if( familycode_j(j)==5 ){
      // familycode = 5 :  Categorical
      if(R_FINITE(asDouble(y_ij(i,j)))){
        if(y_ij(i,j) > 0){
          jnll_ij(i,j) -= log(mu_vj(v_i(i),j));
          sumobs_ig(i,group_j(j)-1) += y_ij(i,j);
        }
      }else{
        sumobs_ig(i,group_j(j)-1) = NAN;
      }
    }
  }}
  jnll += jnll_ij.sum();

  // Likelihood for base level per group
  for(int i=0; i<n_i; i++){
  for(int g=0; g<n_g; g++){
    if(R_FINITE(asDouble(sumobs_ig(i,g)))){
      if( sumobs_ig(i,g) == 0 ){
        jnll_ig(i,g) = -1.0 * log( 1.0 / (1.0 + sumpred_vg(v_i(i),g)));
      }
    }
  }}
  jnll += jnll_ig.sum();

  if( experiments_type == "BH" ){
    DATA_INTEGER( j_logMASPS );
    DATA_IVECTOR( v_p );
    DATA_IVECTOR( p_k );
    DATA_MATRIX( W_kz ); // W_kz.col(0)=R;   W_kz.col(1)=S
    DATA_MATRIX( U_pz ); // U_pz.col(0)=log_SPR0;   U_pz.col(1)=log_M
    PARAMETER_VECTOR( logsigmaR_p );
    PARAMETER_VECTOR( logb_p );
    int n_k = W_kz.rows();
    int n_p = U_pz.rows();

    vector<Type> MLSPS_p( n_p );
    vector<Type> loga_p( n_p );
    for( int p=0; p < n_p; p++ ){
      // Using 1 + value to ensure MLSPS > 1, matching FishLife v2
      MLSPS_p(p) = 1.0 + exp(x_vj( v_p(p)-1, j_logMASPS-1 )) / (1 - exp(-1 * exp(U_pz(p,1))));
      loga_p(p) = log(MLSPS_p(p)) - U_pz(p,0);
    }

    vector<Type> jnll_k( n_k );
    vector<Type> logmu_k( n_k );
    vector<Type> sigmaR_p = exp(logsigmaR_p);
    vector<Type> b_p = exp(logb_p);
    for( int k=0; k < n_k; k++ ){
      logmu_k(k) = loga_p(p_k(k)-1) + log( W_kz(k,1) / (1.0 + W_kz(k,1) / b_p(p_k(k)-1) ) );
      jnll_k(k) = -1 * dnorm( log(W_kz(k,0)), logmu_k(k), sigmaR_p(p_k(k)-1), true );
    }
    jnll += jnll_k.sum();
    REPORT( logmu_k );
    REPORT( MLSPS_p );
    REPORT( jnll_k );

    // Calculate CV in logmu_k by population
    vector<Type> sum_logmu_p(n_p);
    vector<Type> num_logmu_p(n_p);
    vector<Type> var_logmu_p(n_p);
    sum_logmu_p.setZero();
    num_logmu_p.setZero();
    var_logmu_p.setZero();
    // Calculate mean of logmu_k by population
    for( int k=0; k<n_k; k++ ){
      num_logmu_p(p_k(k)-1) += 1;
      sum_logmu_p(p_k(k)-1) += logmu_k(k);
    }
    vector<Type> mean_logmu_p = sum_logmu_p / num_logmu_p;
    // Calculate SD of logmu_k by population
    for( int k=0; k<n_k; k++ ){
      var_logmu_p(p_k(k)-1) += pow( logmu_k(k) - mean_logmu_p(p_k(k)-1), 2 );
    }
    vector<Type> sd_logmu_p = pow( var_logmu_p, 0.5 );

    // Penalize low variance in predictive recruitment
    Type Pen_lowvar_lnRhat = 1;    // Value used in Thorson-2020
    //if( Pen_lowvar_lnRhat > 0 ){
      jnll -= Pen_lowvar_lnRhat * sum(log(sd_logmu_p));
    //}
    REPORT( sd_logmu_p );
  }

  // Calculate intercept
  vector<Type> root_j = x_vj.row(vroot);
  vector<Type> intercept_j( n_j );
  intercept_j = (I_jj - Rho_jj) * root_j.matrix();
  REPORT( root_j );
  REPORT( intercept_j );
  ADREPORT( intercept_j );

  // Extract elements of sparse matrix
  //vector<Type> nonzeroRho_z( Rho_jj.nonZeros() );
  //for( int z = 0; z<Rho_jj.nonZeros(); z++ ){
  //  nonzeroRho_z(z) = Rho_jj.coeffRef(z);
  //}
  vector<Type> nonzeroRho_z = get_nonzero_elements( Rho_jj );

  // Reporting
  REPORT( rho_v );
  REPORT( var_v );
  //REPORT( V_jj );
  REPORT( Q_jj );
  REPORT( Rho_jj );
  REPORT( Gamma_jj );
  REPORT( jnll );
  REPORT( jnll_v );
  REPORT( jnll_ij );
  if( n_g > 0 ){
    REPORT( jnll_ig );
    REPORT( sumobs_ig );
    REPORT( sumpred_vg );
  }
  REPORT( alpha );
  REPORT( x_vj );
  REPORT( yhat_ij );  // Testing for cAIC
  REPORT( mu_vj );  // response-scale predictor (including categorical traits)
  REPORT( eps_vj );
//  ADREPORT( Rho_jj );
  ADREPORT( nonzeroRho_z );
  ADREPORT( alpha );
  ADREPORT( lambda );
  ADREPORT( kappa );
  return jnll;
}
