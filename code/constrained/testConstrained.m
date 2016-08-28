% (c) Luca Nenna
%grid
N=200;
xa=0;
xb=1;
ya=0;
yb=1;
x=linspace(xa,xb,N);
y=linspace(ya,yb,N);

%marginals
mu=@(x) (x>=xa).*(x<=xb);
nu=@(x) (x>=ya).*(x<=yb);
oo=ones(1,N);
lambda=1.e-3;
%cost matrix
h=@(x,y) exp(-(0.5*abs(x-y).^2)/lambda);
gamma0=h(oo'*x,y'*oo);
%capacity constrained
hbar=1.5*ones(N,N);
%if the marginals have not the same total mass they must be normalized
m=mu(x)';
n=nu(y);
[gamma] = ConstrainedOT(x,y,m,n,gamma0,hbar,lambda);
