function [x1,x2,y1,y2,xI,yI,xf,yf,mmP,m1,n1,totalTime,it] = partialOT(tol,epsilon,N,mass,eta)
%PARTIALOT 
%Input
% tol: tolerance for th relative error
% epsilon: parameter for the entropic regularisation
% N: number of grid points
% mass: fraction of the min[mass(mu),mass(nu)] to transport
% eta: threshold to detect the active region
%Output
% x1,x2: coordinate of mu
% y1,y2: coordinate of nu
% xI,yI: points of the Active Source
% xf,yf: points of the Active Target
% m1: first marginal of gamma
% n1: second marginal of gamma
% mmP: matix for the .png
% totalTime: time (s) to converge
% it: number of iterations of Dykstra
% Luca Nenna 12/12/2014

%define the grid
xa=0;
xb=1;
ya=-0.25;
yb=0.55;
dx=0.75;
dy=0.75;
x1=linspace(xa,xb,N);
x2=linspace(xa,xb,N);

y1=linspace(xa,xb,N);
y2=linspace(xa,xb,N);

wx1=(xb-xa)/N;
wy1=(yb-ya)/N;
wx=x1(2)-x1(1);
wy=wy1;
[X1,X2]=meshgrid(x1,x2);
[Y1,Y2]=meshgrid(y1,y2);
%define the density (see density.m)
[mu,nu,areaX,areaY]=density(dx,dy);

[I1,J1,l1]=find(abs(X1)~=inf);
[I1,J1,l2]=find(abs(X2)~=inf);
[I2,J2,g1]=find(abs(Y1)~=inf);
[I2,J2,g2]=find(abs(Y2)~=inf);
x1=X1(l1);
x2=X2(l2);
y1=Y1(g1);
y2=Y2(g2);
[M,M1]=size(X1);
m=double(mu(x1,x2));
n=double(nu(y1,y2));
%we consider only the gridpoints where mu!=0 and nu!=0
[I]=find(m);

[J]=find(n);
m=m(I);
x1=x1(I);
x2=x2(I);
n=n(J);
y1=y1(J);
y2=y2(J);
N1=length(x1);
N2=length(y1);
wx1=areaX/N1;
wy1=areaY/N2;
%quadratic cost
h=@(x1,x2,y1,y2) exp(-((0.5*(abs(x1-y1).^2+abs(x2-y2).^2)))/epsilon);


%fixed mass
mass=mass*min([areaX areaY]);


err=1.0;
count=1;

gamma0=zeros(N1,N2);
%initial point for the Dykstra's Algo
for i=1:N2
gamma0(:,i)=(h(x1,x2,y1(i),y2(i)));
end
g=gamma0;
q1=ones(N1,N2);
q2=q1;
q3=q1;
gamma1=gamma0;
gamma2=gamma0;
gamma3=gamma0;
ga=gamma0;
gammaN=gamma0;
qtmp=q1;
count=1;
C=3;
%% loop
    tic;
while (err)>tol

if mod(count,C)==1
    str = sprintf('Proj on C1');
    disp(str);

    tmp1=exp(log(ga)+log(q1));%gamma.*q1;

    gammaN=exp(log(tmp1.*mass)-log(sum(sum(exp(log((wx1)*tmp1*(wy1)))))));
    qtmp=q1;
    q1=exp(log(tmp1)-log(gammaN));

    gamma1=gammaN;   
 
elseif mod(count,C)==2
    str = sprintf('Proj on C2');
    disp(str);
  
    tmp1=ga.*q2;
    tmp1(isnan(tmp1))=0;
    gammaN=tmp1;
    k=find((sum(tmp1,2)*wy1)>m);

    gammaN(k,:)=tmp1(k,:).*((m(k)./(sum(tmp1(k,:),2)*wy1))*ones(1,N2));

    q2=(tmp1./gammaN);

    gamma2=gammaN;

else
    str = sprintf('Proj on C3');
    disp(str);
    tmp1=ga.*q3; tmp1(isnan(tmp1))=0;
    gammaN=tmp1;
    k=find((sum(tmp1,1)*wx1)'>n);
    gammaN(:,k)=tmp1(:,k).*(ones(N1,1)*(n(k)'./(sum(tmp1(:,k),1)*wy1)));

    
    q3=(tmp1./gammaN);

    gamma3=gammaN;
end
    

if count >=1

  err=(sum(sum(exp(log((wx1)*abs(gammaN-gamma1)*(wy1)))))+sum(sum(exp(log((wx1)*abs(gammaN-gamma2)*(wy1)))))+sum(sum(exp(log((wx1)*abs(gammaN-gamma3)*(wy1))))));%*wy1*wx1 ;
end

    str = sprintf('Error at step (%d): %d', count, err);
    disp(str);
ga=gammaN;
count=count+1;
end
totalTime=toc;
it=count-1;
%% post-process
%compute the marginals of ga
m1=m;
for i=1:N1
    m1(i)=sum(ga(i,:).*wy1');
end
n1=n;
for i=1:N2
    n1(i)=sum(ga(:,i).*wx1);
end


xf=y1(abs(n1)/(mass)>eta);
yf=y2(abs(n1)/(mass)>eta);
xI=x1(abs(m1)/(mass)>eta);%>1.e-4);
yI=x2(abs(m1)/(mass)>eta);%>1.e-4);

mm=double(mu(X1,X2));
mmP=zeros(size(mm));

for i=1:length(x1)
j=find(X1==x1(i) & X2==x2(i));
mmP(j)=1;
end

for i=1:length(y1)
j=find(Y1==y1(i) & Y2==y2(i));
mmP(j)=1;
end

for i=1:length(xI)
j=find(X1==xI(i) & X2==yI(i));
mmP(j)=3;
end

for i=1:length(xf)
j=find(Y1==xf(i) & Y2==yf(i));
mmP(j)=2;
end



end

