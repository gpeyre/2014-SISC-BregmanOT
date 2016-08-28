function [c,MVP,objectives] = bregmanWassersteinBarycenter(C,M,iterations,lambda,useGPU,toleranceDifference,weights)
% INPUT:
% C = d x N , N histograms of size d
% M = ground metric.
% iterations = total number of iterations
%
% ADVICE: divide M by median(M) to have a natural scale
% for lambda
objectives=[];

n=size(C,2);
d=size(M,1);

if nargin<3,
    iterations=100;
end

if nargin<5,
    useGPU=false;
end

if nargin<6,
    toleranceDifference=1e-3;
end

if nargin<7,
    weights=ones(n,1)/n;
else
    weights=reshape(weights,n,1); % just to make sure we get a column vector.
end



K=exp(-lambda*M);

if useGPU    
    C=gpuArray(C);
    K=gpuArray(K);    
end

K(K<1e-300)=1e-300;

compt=0;
differ=inf;
objectives=[];
matrixVector=0;
MVP=[];


% two first projections are simpler because u is necessarily a matrix of ones,
% simplifying computations.
UKv=K*(bsxfun(@rdivide,C,(sum(K))')); % first iteration, U = ones.
u=bsxfun(@ldivide,UKv,exp(log(UKv)*weights));

% loop below 
while compt<iterations && differ>toleranceDifference,
    % projection sur les N contraintes de marge fixees
    % le U ci-dessous represente un D(u), i.e. 
    % UKv represente la matrice des D(u_k) K v_k.        
    UKv=u.*(K*(C./(K'*u)));
    matrixVector=matrixVector+2;
    %MVP=[MVP,matrixVector];
    compt=compt+1;    
    % projection sur l'egalite des N marges variables    
    %u=u.*repmat(exp(mean(log(UKv),2)),1,size(C,2))./UKv;      
    u=bsxfun(@times,u,exp(log(UKv)*weights))./UKv;          
    matrixVector=matrixVector+1; % summing operation is quadratic. 
    % log is also quite time-consuming. exp is only carried out on a vector
    MVP=[MVP,matrixVector];
    
    
    % what follows has only a marginal use and can be commented out.
    differ=sum(std(UKv,1,2))
    objectives=[objectives,differ];
    
    if mod(compt,5)==1,
        differ=sum(std(UKv,1,2))    
        compt        
        c=mean(UKv,2);
%         subaxis(6,10,6,2,5,5,'Spacing', 0.005, 'Padding', 0, 'Margin', 0.005)
%         imagesc(reshape(c,round(sqrt(length(c))),round(sqrt(length(c)))));        
%         %imagesc(reshape(c,63,64));        
%         box on;
%         set(gca,'xtick',[],'ytick',[])
%         xlabel('');
%         ylabel('');
%         title(num2str(compt));
%         pause(.01);
    end
end

c=mean(UKv,2);



