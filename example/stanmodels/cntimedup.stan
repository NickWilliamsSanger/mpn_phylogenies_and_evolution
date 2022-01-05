data{
  int N; //number of dup muts (vaf=2/3)
  int M; //number of non dup muts (vaf=1/3)
  real pdup; //sensitivity of capturing dup
  real pndup; //sensitivity of capturing ndup
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
  N ~ poisson(l*x*0.5*pdup);
  M ~ poisson((l*x*0.5 + l*1.5*(1-x))*pndup);
}

generated quantities {
  real t;
  t=l*x;
}
