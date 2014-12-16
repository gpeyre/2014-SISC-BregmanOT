function [gamma] = ConstrainedOT(x,y,m,n,gamma0,hbar,lambda)
%CONSTRAINEDOT This function solves the capacity constrained transport by using the Dykstra Algorithm
%see section 8.2 of  Benamou,Carlier,Cuturi,Nenna and Peyré,  ¨Iterative Bregman Projections for Regularized Transportation Problems¨ 
%INPUT
%   x,y grid associated to the first and the second marginal
%   m,n are \mu and \nu
%   gamma0 must be the matrix \gamma_{ij}=exp^{-c_{ij}/lambda} where c_{ij} is the cost matrix for the OT problem
%   hbar is the matrix associated to the capacity constraint on \gamma 
%   lambda is the entropic regularization
%OUTPUT
% gamma is the optimal transport plan

% (c) Luca Nenna



N=length(x);
wx=x(2)-x(1);
wy=y(2)-y(1);

err=1.0;



%matrices for Dykstra
q1=ones(N,N);
q2=q1;
q3=q1;
gamma1=gamma0;
gamma2=gamma0;
gamma3=gamma0;
gamma=gamma0;
gammaN=gamma;
oo=ones(1,N);
count=1;

C=3;
%% loop

while (err)>1.e-4 

if mod(count,C)==1
    str = sprintf('Proj on C1');
    disp(str);
    tmp1=gamma.*q1;

    gammaN=min(tmp1,hbar);
    q1=(tmp1./gammaN);
    gamma1=gammaN;
elseif mod(count,C)==2
    str = sprintf('Proj on C2');
    disp(str);
    tmp1=gamma.*q2;
    gammaN=tmp1.*((m./(sum(wy*tmp1,2)))*oo);
    q2=(tmp1./gammaN);
    gamma2=gammaN;

else
    str = sprintf('Proj on C3');
    disp(str);
    tmp1=gamma.*q3;
    gammaN=tmp1.*(oo'*(n./(sum(wx*tmp1,1))));    
    q3=(tmp1./gammaN);
    gamma3=gammaN;
end
    

if count >=3
err=(sum(sum(abs(gammaN-gamma1)))+sum(sum(abs(gammaN-gamma2)))+sum(sum(abs(gammaN-gamma3))))*wx*wy;
end

str = sprintf('Error at step (%d): %d', count, err);
disp(str);
gamma=gammaN;
count=count+1;
end


%plot gamma
pcolor(x,y,gamma)
shading flat;
colormap gray
colormap(flipud(colormap))
colorbar;


end

