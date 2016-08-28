function [Call,MVP,objectives] = smoothConjugateDiscreteWassersteinBarycenter(C,M,iterations,t0,lambda,useGPU,toleranceDifference,StartingTrick)
% INPUT:
% R = d x N , N histograms of size d
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

if nargin<4,
    t0=10;
end

if nargin<6,
    useGPU=false;
end

if nargin<7,
    toleranceDifference=1e-4;
end


K=exp(-lambda*M);

if useGPU
    M=gpuArray(M);
    C=gpuArray(C);
    K=gpuArray(K);    
end
 
%StartingTrick=0;

totalMatrixVectorProducts=0;

K(K<1e-200)=1e-200;
if ~StartingTrick, % initialisation naive
    G=zeros(d,n);
    U=C; % first iteration is trivial in this case, U= database itself.
    G=G+t0*U/sqrt(1);
    %figure(3); imagesc(G);
    G=G-repmat(mean(G,2),1,n);
else
    if useGPU,
        unif=gpuArray(ones(d,1)/d);
    else
        %unif=ones(d,1)/d;
        unif=rand(d,1);
        unif=unif/sum(unif);
    end
    [~,u,~,matrixVectorProducts]=sinkhornTransportHotStartG(unif,C,K,K.*M,[],'marginalDifference',inf,.4,useGPU);
    totalMatrixVectorProducts=totalMatrixVectorProducts+matrixVectorProducts;
    %[~,~,u,~]=sinkhornTransport(unif,C,K,K.*M,lambda,[],[],[],[],[],[],[0,0,1,0]);
    G=log(u)/lambda;
    %    G=G-repmat(mean(G,1),n,1); % center gradient.
    G=bsxfun(@minus,G,mean(G,2)); % center gradient.
end

if useGPU
    G=gpuArray(G);
end

compt=1;

cold=ones(d,1);
c=zeros(d,1);
differ=inf;

Call=[];
MVP=[];
objectives=[];

LAG=4;
while compt<iterations && differ>toleranceDifference
    %lambda=lambda+100/iterations
    
    %pause();
    expG=exp(lambda*G);
    U=expG.*( K* (C./(K*expG))); % compute gradient
    totalMatrixVectorProducts=totalMatrixVectorProducts+2;
    %     sum(G.*U)+... % compute objective <g,u>
    %     sum(sum(expG.*(M.*(K*(C./(K*expG)))))+...% transportation costs
    %
    %disp(['Sum of conjugate objectives:',num2str(sum(-Dist))]);
    
    G=G-t0*U;%/sqrt(compt); % remove gradient
    G=bsxfun(@minus,G,mean(G,2)); % center gradient.
    %     expG=exp(lambda*G);
    %
    %     [D,~,u,~]=sinkhornTransport(U,C,K,K.*M,lambda,[],[],[],[],[],[],[1,0,1,0],expG);
    %     disp(['Sum of primal objective:',num2str(sum(D))]);
    %     G=log(u)/lambda;
    %     G=G-repmat(mean(G,1),d,1); % center gradient.
    %     G=G-repmat(mean(G,2),1,n); % center gradient.
    %figure(3); imagesc(G); title('G');
    %figure(4); imagesc(U); title('U');
    
    %cold=c;
    c=mean(U,2);
    differ=norm(c-cold,1);
      
    
    if mod(compt,LAG)==1 || compt==iterations-1,
        LAG=2*LAG;
        Call=[Call,c];
        MVP=[MVP,totalMatrixVectorProducts];
        objectives=[objectives,sum(std(U,0,2))];
    end
    
    
%     if 1, %compt<250
%         figure(1);plot(c); title(num2str(compt)); 
%         pause(.01)
%     else
%         pause();
%     end
    
    if 0
        figure(2);
        imagesc(reshape(mean(U,2),round(sqrt(d)),round(sqrt(d))));
        colorbar;
        box on;
        set(gca,'xtick',[],'ytick',[])
        xlabel('');
        ylabel('');
    end
    
    if 0 || mod(compt,1)==1,
        compt
        %subaxis(6,10,6,2,5,5,'Spacing', 0.005, 'Padding', 0, 'Margin', 0.005)
        figure(3);
        imagesc(reshape(c,round(sqrt(length(c))),round(sqrt(length(c)))));        
        box on;
        set(gca,'xtick',[],'ytick',[])
        xlabel('');
        ylabel('');
        title(num2str(compt));
        pause(.1);
    end
    compt=compt+1;
end




