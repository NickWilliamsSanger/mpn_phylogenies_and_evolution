data{
  int N; //number of homs
  int M; //number of hets
  real phet; //probability of observing het
  real phom; //probability of observing hom
  real L; //Expected number of mutations
  real Lstd; //Standard Error in L
}

parameters {
  real<lower=0,upper=1> x;
  real<lower=1> l;
}

model {
  x ~ uniform(0,1);
  l ~ normal(L,Lstd);
  N ~ poisson(l*x*phom);
  M ~ poisson(l*(1-x)*phet);
}

generated quantities {
  real t;
  t=l*x;
}
