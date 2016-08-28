
N=200 ; 
Nmarg= 4^2; %use 16 marginals otherwise change the function computeGamma
map=3;
lambda=5.e-4;

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
    mu=m(x)'/N;


gamma0=zeros(N,N); %this is the matrix c^{i_{j}i_{j+1}} j=1,...,Nmarg-1
gamma1=zeros(N,N); %this is the matrix c^{Nmarg1}
for i=1:N


gamma0(i,:)=h(x(i),y);
gamma1(i,:)=h(t,x(i));

oo = ones(N,1) ;

end
gamma0=sparse(gamma0);
gamma1=sparse(gamma1);
A=ones(N,Nmarg);
% A(:,1)=mu;  a quoi ca sert 
err=1;

count=1;
while err>1.e-4

A_old = A ;


%% optimize - precompute gamma_0^N
    
for i=1:Nmarg
B = A ;     
B(:,i) = oo ; % remove one kp
B1 = eye(N,N) ;  
rr = circshift([1:Nmarg]',1-i)'  ; %circshift to construct Gamma start sum at i+1
for j = rr
gamma = gamma0 ; 
if (j==Nmarg) 
    gamma = gamma1 ; 
end 
    B1 = B1*diag(B(:,j))*gamma ; 
end  

A(:,i)=(mu./diag(B1));

end

err=norm((A_old-A),Inf)/norm(A,Inf) ;
str = sprintf('Error at step (%d): %s', count, err);
disp(str);
count=count+1;

end

figure(1) 
l = sqrt(Nmarg) ; 

%% 2marginals visu 
err=zeros(length(2:Nmarg-1),1);
count=1;
for i=2:Nmarg;%-1 ; 
B1 =  eye(N,N) ; 
for j = 1:i-1 
B1 = B1*diag(A(:,j))*gamma0 ; 
end   
B2 = eye(N,N) ;  
for j = i:Nmarg-1  
B2 = B2*diag(A(:,j))*gamma0 ; 
end  
B2 = B2*diag(A(:,Nmarg))*gamma1 ; 

%% product of B1 and B2 component wize 
%r=sum(((B1.*B2')'),2);
gammaN=((B1.*B2')');%/diag(r);
%gammaN=gammaN/diag(sum(gammaN,1));
subplot(l,l,i);
pcolor(x,x,gammaN);

t=(i-1)/Nmarg;
gammaEx=solutionEuler(t,N);
err(count)=norm(gammaN-gammaEx,Inf);
axis square 
shading flat;
caxis([0 0.002])
count=count+1;

end 
colorbar 


