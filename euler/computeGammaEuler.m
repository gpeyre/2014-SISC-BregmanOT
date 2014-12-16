N=200;%N nuber of points 
Np=50;%Np number of particles 
%t time 
v=@(A,omega) pi*sqrt(A.*(1-A)./2).*cos(pi*omega);
G=@(t,A,omega) 0.5+(A-0.5).*cos(pi.*t)+v(A,omega).*sin(pi.*t);
m=@(x) (x>=0).*(x<=1);

l=linspace(0,1,N);

gamma=zeros(N,N);
omega=linspace(1,Np,Np)/Np;

[A,O]=meshgrid(l,omega);


figure
for k=1:9
t=(k-1)/8;
G1=G(t,A,O);

a1=max(max(G1,[],1));
a2=min(min(G1,[],1));
dx=1/N;
len=a1-a2;
for i=1:N
H=(G(t,l(i),omega)-a2);%G>=0
for j=1:N
h=find(H>=((j-1)*dx*len)& (H<j*dx*len))/(Np*N); %normalization \gamma_{i,j}\in[0,1]
gamma(i,j)=length(h);
end
end

subplot(3,3,k)
pcolor(l,l,gamma)
shading flat
title(['T= ',num2str(t)])
colormap gray
colormap(flipud(colormap))
caxis([0 0.0005]);
end
