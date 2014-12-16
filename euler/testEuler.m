N=200 ; 
Nmarg= 4^2; %use 16 marginals otherwise change the function computeGamma
map=3;
lambda=5.e-4;


[x,phi,err] = IPFPEuler(lambda,N,Nmarg,map);
