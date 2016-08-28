function [D,l,m,matrixVectorProducts]=sinkhornTransportHotStartG(a,b,K,U,l,stoppingCriterion,p_norm,tolerance,useGPU)
% Compute a vector of N Sinkhorn distances, from a single histogram a
% to N histograms b_1, ... , b_N.  Outputs a vector of distances
% D= [d(a,b_1), d(a,b_2), ... , d(a,b_N)].
% If needed, the function also outputs diagonal scalings to recover smoothed optimal
% transport between each of the pairs (a,b_i).
%
%---------------------------
% Required Inputs:
%---------------------------
% a is a d1 vector in the probability simplex (nonnegative, summing to one)
% b is a d2 x N matrix of N vectors in the probability simplex
% K is a d1 x d2 matrix, equal to exp(-lambda M), where M is the d1 x d2
% matrix of pairwise distances between bins of a and all the bins in the b_1,...b_N histograms.
% U = K.*M is a d1 x d2 matrix, pre-stored for efficiency
%
%
%---------------------------
% Optional Inputs:
%---------------------------
% stoppingCriterion in {'marginalDifference','distanceRelativeDecrease'}
%   - marginalDifference (Default) : checks whether the difference between
%              the marginals of the current optimal transport and the
%              theoretical marginals set by a b_1,...,b_N are satisfied.
%   - distanceRelativeDecrease : only focus on convergence of the vector
%              of distances
%
% p_norm: parameter in {(1,+infty]} used to compute a stoppingCriterion statistic
% from N numbers (these N numbers might be the 1-norm of marginal
% differences or the vector of distances.
%
% tolerance >0
%
%
%---------------------------
% Output
%---------------------------
% D : vector of N dual-sinkhorn divergences
% l : d1 x N matrix of left scalings
% m : d2 x N matrix of right scalings
%
% The smoothed optimal transport between (a,b_i) can be recovered as
% T_i = diag(l(:,i)) * K * diag(m(:,i));
%
% or, equivalently and substantially faster:
% T_i = bsxfun(@times,m(:,i)',(bsxfun(@times,l(:,i),K)))
%
%
% Relevant paper:
% M. Cuturi,
% Sinkhorn Distances : Lightspeed Computation of Optimal Transport,
% Advances in Neural Information Processing Systems (NIPS) 26, 2013

% This code, (c) Marco Cuturi 2013 (see license block below)
% v0.0 first version 20/11/2013

%% Processing optional inputs

if nargin<6,
    stoppingCriterion='marginalDifference'; % check marginal constraints by default
end

if nargin<7,
    p_norm=inf;
end

if nargin<8,
    tolerance=1e-3;
end

if nargin<9,
    useGPU=false;
end




%% Check that a only has >0 components
I=(a>0);

if sum(I)<length(a), % need to update some vectors and matrices if a does not have full support
    K=K(I,:);
    if ~isempty(U), U=U(I,:); end
    a=a(I);
    if ~isempty(l), l=l(I,:); end
end


%% Initialization
maxIter=10000;
compt=0;

matrixVectorProducts=0;
if nargin<5 || isempty(l),
    if useGPU,
        l=ones(length(a),size(b,2),'gpuArray')/length(a);
    else
        l=ones(length(a),size(b,2))/length(a);
    end
    
end
rK=bsxfun(@rdivide,K,a); % precomputation of this matrix saves a d1 x N Schur product
matrixVectorProducts=matrixVectorProducts+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(stoppingCriterion,'distanceRelativeDecrease')
    Dold=ones(1,size(b,2)); %initialization of vector of distances.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checking dimensionality:
if size(b,2)>size(b,1),
    BIGN=true;
else
    BIGN=false;
end



while compt<maxIter,
    
    if BIGN
        l=(1./(rK*(b./(K'*(1./(rK*(b./(K'*(1./(rK*(b./(K'*(1./(rK*(b./(K'*l)))))))))))))))); % main iteration of Sinkhorn's algorithm
        matrixVectorProducts=matrixVectorProducts+8;            
%         l=1./(rK*(b./(K'*l))); % main iteration of Sinkhorn's algorithm
%         matrixVectorProducts=matrixVectorProducts+2;
    else
        l=1./(rK*(b./((1./(rK*(b./((1./(rK*(b./((1./(rK*(b./(l'*K)')))'*K)')))'*K)')))'*K)'));
        matrixVectorProducts=matrixVectorProducts+8;            
%         l=1./(rK*(b./(l'*K)'));
%         matrixVectorProducts=matrixVectorProducts+2;
    end
    
    %l=1./(rK*(b./(K'*l))); % main iteration of Sinkhorn's algorithm
    
    if 1, %mod(compt,1)==1,   % check every 10 Sinkhorn iterations.
        % split computations to recover right and left scalings.
        m=b./(K'*l);
        l=1./(rK*m);
        matrixVectorProducts=matrixVectorProducts+2;
        % check stopping criterion
        switch stoppingCriterion,
            case 'distanceRelativeDecrease',
                D=sum(l.*(U*m));
                if norm(D./Dold-1,p_norm)<tolerance,
                    break;
                end
                Dold=D;
                
            case 'marginalDifference',
                toltest=norm(sum(abs(m.*(K'*l)-b)),p_norm);
                %disp([num2str(compt),' - ',num2str(toltest)]);
                if toltest<tolerance, % norm of all || . ||_1 differences between the marginal of the current solution with the actual marginals.
                    break;
                end
        end
        compt=compt+1;
    end
    compt=compt+1;
end

if compt==maxIter, % if did not converge return NaN and display problem.
    D=nan(1,size(b,2));
    disp('Sinkhorn Transport Reached the Maximal Number of Iterations --- ');
elseif strcmp(stoppingCriterion,'marginalDifference') && ~isempty(U), % if we have been watching marginal differences, we need to compute the vector of distances.
    D=sum(l.*(U*m));
else
    D=nan;
end

if nargout>1, % user wants scalings
    L=l;
    if useGPU,
        l=gpuArray(zeros(length(I),size(b,2)));
    else
        l=zeros(length(I),size(b,2));
    end
    l(I,:)=L;
end

% ***** BEGIN LICENSE BLOCK *****
%  * Version: MPL 1.1/GPL 2.0/LGPL 2.1
%  *
%  * The contents of this file are subject to the Mozilla Public License Version
%  * 1.1 (the "License"); you may not use this file except in compliance with
%  * the License. You may obtain a copy of the License at
%  * http://www.mozilla.org/MPL/
%  *
%  * Software distributed under the License is distributed on an "AS IS" basis,
%  * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
%  * for the specific language governing rights and limitations under the
%  * License.
%  *
%  * The Original Code is Sinkhorn Transport, (C) 2013, Marco Cuturi
%  *
%  * The Initial Developers of the Original Code is
%  *
%  * Marco Cuturi   mcuturi@i.kyoto-u.ac.jp
%  *
%  * Portions created by the Initial Developers are
%  * Copyright (C) 2013 the Initial Developers. All Rights Reserved.
%  *
%  *
%  ***** END LICENSE BLOCK *****
