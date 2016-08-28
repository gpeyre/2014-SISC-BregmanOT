function c = WassersteinBarycenter_lesong(C,M,iterations,lambda,useGPU,toleranceDifference,weights)
% INPUT:
% C = d x N , N histograms of size d
% M = ground metric.
% iterations = total number of iterations
% lambda = entropic regularization
% useGPU = binary flag.
% toleranceDifference = problem is convex. When iterates change by less
% than this value in l_1 norm the algorithm stops.
% weights = this program computes the isobarycenter by default, a weighted
% one if this parameter is provided.


% OUTPUT
% c : barycenter according to weights
% ADVICE: divide M by median(M) to have a natural scale
% for lambda

% (c) 2014 Marco Cuturi, Gabriel Peyré, Guillaume Carlier

n=size(C,2);


if nargin<3,
    iterations=10000;
end

if nargin<4,
    lambda=100;
end


if nargin<5,
    useGPU=false;
end

if nargin<6,
    toleranceDifference=1e-4;
end

if nargin<7,
    weights=false;
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


% two first projections are simpler because u is necessarily a matrix of ones,
% simplifying computations.
UKv=K*(bsxfun(@rdivide,C,(sum(K))')); % first iteration, U = ones.
u=bsxfun(@ldivide,UKv,exp(mean(log(UKv),2)));

% loop below 
while compt<iterations && differ>toleranceDifference,
    compt=compt+1;    
    % projection sur les N contraintes de marge fixees
    % le U ci-dessous represente un D(u), i.e. 
    % UKv represente la matrice des D(u_k) K v_k.        
    UKv=u.*(K*(C./(K'*u)));
    
    % projection sur l'egalite des N marges variables       
    if ~weights, % isobarycenter
        u=bsxfun(@times,u,exp(mean(log(UKv),2)))./UKv;      
    else % weights
        u=bsxfun(@times,u,exp(log(UKv)*weights))./UKv;      
    end
    
    
    % what follows has only a marginal use and can be commented out.        
    if mod(compt,10)==1,
        differ=sum(std(UKv,1,2));
        disp(['Iteration ',num2str(compt),' - ',num2str(differ)]);
        % put anything here if you want to plot current solution
    end
end

c=mean(UKv,2);



