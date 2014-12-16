function [x,phi,err] = IPFPEuler(lambda,N,Nmarg,map)
%USE ONLY 16 MARGINALS otherwise change function computeGamma
%lambda is the entropic regularization
% Nmarg is the number of marginals
%N is the number of points for the grid
% map is a measure preserving map 
%     1:S(x)=min(2*x,2-2*x);
%     2:S(x)=mod(x+0.5,1);
%     3:S(x)=1-x;
% (c) Luca Nenna

x=linspace(0,1.,N);
y=x;


    S1  = @(x) min(2*x,2-2*x);
    S2  = @(x) mod(x+0.5,1);
    S3  = @(x) 1-x;
    c  = @(x,y) 0.5*abs(x-y).^2;

    h=@(x,y) exp(-c(x,y)/lambda);
    m=@(x) (x>=0).*(x<=1);

    if map==1 
    t=S1(x);
    elseif map==2
    t=S2(x);
    else
    t= S3(x);
    end
    mu=m(x)';




gamma0=zeros(N,N); %this is the matrix c^{i_{j}i_{j+1}} j=1,...,Nmarg-1
gamma1=zeros(N,N); %this is the matrix c^{Nmarg1}
for i=1:N


gamma0(i,:)=h(x(i),y);
gamma1(i,:)=h(t,x(i));


end

A=ones(N,Nmarg);
err=1;

count=1;
while err>1.e-5 

A(:,end)=projection(mu,A(:,1:end-1),gamma0,gamma1,Nmarg,Nmarg);
for i=1:Nmarg-1
A(:,i)=projection(mu,[A(:,1:i-1) A(:,i+1:end)],gamma0,gamma1,i,Nmarg);
end

an=projection(mu,A(:,2:end),gamma0,gamma1,1,Nmarg);

err=norm((an-A(:,1)),Inf)/norm(abs(A(:,1)),Inf);
str = sprintf('Error at step (%d): %s', count, err);
disp(str);
A(:,1)=an;
count=count+1;
end

phi=(log(A(:,1)))*lambda;

[n,err]=computeGamma(A,gamma0,gamma1,x,map);


end



function [b]=projection(mu,B,gamma0,gamma1,pot,Nmarg)
[m n]=size(B);
[m1 n1]=size(gamma0);
B1=eye(m1,n1);
%sum over i_{j}!=i_{1}
if pot==1
B1=gamma0;
for i=1:n-1
B1=B1*(diag(B(:,i))*gamma0);
end
%B1=gamma0'.*B1;
B1=B1*(diag(B(:,end))*gamma1);
%sum over i_{j}!=i_{pot}
elseif pot~=1 && pot~=Nmarg

%B1=gamma0;%(diag(B(:,pot-1))*gamma0);

for i=1:n

if n-i+1==(Nmarg-pot)
B1=(diag(B(:,end))*gamma1)*B1;
%n-i+1
elseif  (n-i+1)>(Nmarg-pot)
%pot-i
B1=(diag(B(:,pot-i))*gamma0)*B1;

elseif  (n-i+1)<(Nmarg-pot)

B1=(diag(B(:,n-i+pot))*gamma0)*B1;
%n-i+pot
end

end
%B1=gamma0'.*B1;
B1=gamma0*B1;
%sum over i_{j}!=i_{Nmarg}
else
B1=gamma1;
for i=1:n-1
B1=B1*(diag(B(:,i))*gamma0);
end
%B1=gamma1'.*B1;
B1=B1*diag(B(:,end))*(gamma0);
end

%b=(mu./sum(B1,1)');
b=(mu./diag(B1));
%b(isnan(b))=1.e-300;
%b(isinf(b))=1.e300;

end

function [b]=gammaProj(C,B,D,gamma0,gamma1,pot,x,map)
[m n]=size(B);
[m2 n2]=size(D);
[m1 n1]=size(gamma0);
b=zeros(m1,n1);

%B1=eye(m1,n1);
if pot==0
A=C(:,1);
for i=1:m1
B1=diag(B(:,end))*gamma1;
B1=B1(:,i);
j=n-1;
while j>=1
B1=(diag(B(:,j))*gamma0)*B1;
j=j-1;
end
B1=((C(:,2)*ones(1,m)).*gamma0)*B1;
%B1=B1.*gamma1';
s=gamma0(i,:).*B1';
b(i,:)=A(i,:).*s;
end
%b=((C(:,1)*ones(1,m))'.*(C(:,2)*ones(1,m))'.*gamma0').*B1;



elseif pot==1
A=C(:,2);
for i=1:m1
B1=diag(B(:,end))*gamma0;
B1=B1(:,i);
j=n-1;
while j>=1
B1=(diag(B(:,j))*gamma0)*B1;
j=j-1;
end
B1=(diag(C(:,1))*gamma0)*B1;
%B1=B1.*gamma1';
s=gamma1(i,:).*B1';
b(i,:)=A(i,:).*s;
end



%compute a sum over i_{j}!=i_{1}
elseif pot==4


B1=(diag(B(:,end))*gamma1);
for i=1:n-1
B1=(diag(B(:,n-i))*gamma0)*B1;
end
B1=gamma0*B1;
%s=sum(B1,1);
s=diag(B1);
for i=1:m1
b(i,i)=C(i)*s(i);
end





%compute the transport plan at the final T \gamma_{ij} where j=S(i) (S is the given map)
elseif pot==5
if map==1
S  = @(x) min(2*x,2*m1-2*x+1);
elseif map==2
S = @(x) mod(x+0.5*m1,m1)+1;
end





B1=(diag(B(:,end))*gamma0);
for i=1:n-1
B1=(diag(B(:,n-i))*gamma0)*B1;
end
B1=gamma1*B1;
s=diag(B1);

if map==1 || map==2
for i=1:m1

j=S(i);
b(i,j)=(C(i)*s(i));

end
else
for i=1:m1

b(i,m1-i+1)=(C(i)*s(i));
%b(i,i)=S(x(i));
end
end


else


for k=1:n1
E1=(diag(B(:,end))*gamma0);
E2=(diag(C(:,2))*gamma0); %matrix associated to the j potential in C(:,2)
E=E1(:,k)*E2(k,:);        %we take the k-th component of x_{j}

B1=diag(D(:,end))*gamma1;

j=n2-1;
while j>=1
B1=(diag(D(:,j))*gamma0)*B1;  %sum over all indices i_{z}>i_{j}
j=j-1;
end
E=E*B1;
B1=E;
clear E;
for j=1:n-1                   %sum over all indices i_{z}<i_{j}  
B1=(diag(B(:,n-j))*gamma0)*B1;
end

s=(diag(C(:,1))*gamma0)*B1;

b(:,k)=diag(s);
end



end







end





%here we compute the transport plan at time T=j-1/16


function [n,err]=computeGamma(A,gamma0,gamma1,x,map)


[m n]=size(A);
Nmarg=n;
N=m;
err=zeros(8,1);
%-T-1------------------------------------------------------------------------------
figure
[gamma]=gammaProj([A(:,1)],A(:,2:end),[],gamma0,gamma1,4);

subplot(3,3,1);
pcolor(x,x,gamma);
shading flat;

caxis([0 0.05]);
title(['T=0']);
%legend;
t=0;
gammaEx=solutionEuler(t,N);
err(1)=norm(gamma-gammaEx,Inf);





%-T-1/8------------------------------------

[gamma]=gammaProj([A(:,1) A(:,3)],[A(:,2)],[A(:,4:end)],gamma0,gamma1,2);

subplot(3,3,2);
pcolor(x,x,gamma');
shading flat;

caxis([0 0.05]);
title(['T=1/8']);
t=(3-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(2)=norm(gamma-gammaEx,Inf);
%legend;
%-T-1/4------------------------------------

[gamma]=gammaProj([A(:,1) A(:,5)],[A(:,2:4)],[ A(:,6:end)],gamma0,gamma1,2);

subplot(3,3,3);
pcolor(x,x,gamma');

shading flat;

caxis([0 0.05]);
title(['T=1/4']);
t=(5-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(3)=norm(gamma-gammaEx,Inf);
%legend;
%-T-3/8------------------------------------

[gamma]=gammaProj([A(:,1) A(:,7)],[A(:,2:6)],[ A(:,8:end)],gamma0,gamma1,2);

subplot(3,3,4);
pcolor(x,x,gamma');

shading flat;

caxis([0 0.05]);
title(['T=3/8']);
t=(7-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(4)=norm(gamma-gammaEx,Inf);
%legend;
%-T-1/2------------------------------------

[gamma]=gammaProj([A(:,1) A(:,9)],[A(:,2:8)],[ A(:,10:end)],gamma0,gamma1,2);

subplot(3,3,5);
pcolor(x,x,gamma');

shading flat;

caxis([0 0.05]);
title(['T=1/2']);
%legend;
%-T-5/8------------------------------------
[gamma]=gammaProj([A(:,1) A(:,11)],[A(:,2:10)],[ A(:,12:end)],gamma0,gamma1,2);

subplot(3,3,6);
pcolor(x,x,gamma'); 
shading flat;

caxis([0 0.05]);
title(['T=5/8']);
t=(11-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(5)=norm(gamma-gammaEx,Inf);
%legend;
%-T-3/4------------------------------------
[gamma]=gammaProj([A(:,1) A(:,13)],[A(:,2:12)],[ A(:,14:end)],gamma0,gamma1,2);

subplot(3,3,7);
pcolor(x,x,gamma');  
shading flat;
caxis([0 0.05]);


title(['T=3/4']);
t=(13-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(6)=norm(gamma-gammaEx,Inf);
%legend;
%-T-7/8------------------------------------
[gamma]=gammaProj([A(:,1) A(:,15)],A(:,2:14),A(:,end),gamma0,gamma1,2);
subplot(3,3,8);
pcolor(x,x,gamma');  
shading flat;

caxis([0 0.05]);
title(['T=7/8']);
t=(15-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(7)=norm(gamma-gammaEx,Inf);
%legend;
%-T-8------------------------------------
[gamma]=gammaProj([A(:,16)],A(:,1:end-1),[],gamma0,gamma1,5,x,map);

subplot(3,3,9);
pcolor(x,x,gamma');
shading flat;

caxis([0 0.05]);
title(['T=1']);
t=1;
gammaEx=solutionEuler(t,N);
err(8)=norm(gamma-gammaEx,Inf);
%--------------------------------------
colormap gray
colormap(flipud(colormap))

end

