functions {
  matrix fill_matrix(real absent_nonbasic, real nonbasic_absent, real noIDCC_IDCC, real IDCC_noIDCC, real nonbasic_basic, real basic_nonbasic, int J) {
    matrix[J,J] Q = rep_matrix(0,J,J);
    Q[1,2] = absent_nonbasic;
    Q[2,1] = nonbasic_absent;
    Q[1,3] = absent_nonbasic;
    Q[3,1] = nonbasic_absent;
    Q[2,3] = noIDCC_IDCC;
    Q[3,2] = IDCC_noIDCC;
    Q[2,4] = nonbasic_basic;
    Q[4,2] = basic_nonbasic;
    Q[3,5] = nonbasic_basic;
    Q[5,3] = basic_nonbasic;
    Q[4,5] = noIDCC_IDCC;
    Q[5,4] = IDCC_noIDCC;
    for (j in 1:J) {
      Q[j,j] = -sum(Q[j,]);
    }
    return(Q);
  }
  //compute likelihood via Felsenstein's Pruning Algorithm
  real pruning_vec(int N, int B, int J, int[] child, int[] parent, real[] brlen, matrix tiplik, matrix Q) {
    matrix[N,J] lambda;                  //likelihoods at tips+nodes
    vector[J] pi;                         //stationary probability
    lambda = log(tiplik);
    for (b in 1:B) {
      matrix[J,J] P = matrix_exp(Q*brlen[b]); //via matrix exponentiation
      for (d in 1:J) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    pi = -log(J) + to_vector(lambda[parent[B],]);
    return(log_sum_exp(pi));
  }
}
data {
  int<lower=1> N;                          //number of nodes in large families
  int<lower=1> B;                          //number of branches
  int<lower=1> J;
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];               //length of each branch
  matrix[N,J] tiplik;                 //likelihoods for data at tips+internal nodes in of each tree
}
parameters {
  real<lower=0> absent_nonbasic;
  real<lower=0> nonbasic_absent;
  real<lower=0> noIDCC_IDCC;
  real<lower=0> IDCC_noIDCC;
  real<lower=0> nonbasic_basic;
  real<lower=0> basic_nonbasic;
}
transformed parameters {
  matrix[J,J] Q = fill_matrix(absent_nonbasic, nonbasic_absent, noIDCC_IDCC, IDCC_noIDCC, nonbasic_basic, basic_nonbasic, J);
  real log_lik = pruning_vec(N,B,J,child,parent,brlen,tiplik,Q);
}
model {
  absent_nonbasic ~ gamma(1,1);
  nonbasic_absent ~ gamma(1,1);
  noIDCC_IDCC ~ gamma(1,1);
  IDCC_noIDCC ~ gamma(1,1);
  nonbasic_basic ~ gamma(1,1);
  basic_nonbasic ~ gamma(1,1);
  target += log_lik;
}
