function r = compute_radon_transform(f,theta,transposed)

% compute_radon_transform - compute a fully discrete Radon transform
%
%   Compute either forward transform from f of size (N,N) to (N,Q) arrays, 
%   where Q=length(theta):
%       r = compute_radon_transform(f,theta,0);
%   or transposed (backprojection) 
%       f = compute_radon_transform(r,theta,1);
%
%   Copyright (c) 2014 Gabriel Peyre

N = size(f,1);
Q = length(theta);

if Q>1
    if transposed==0
        r = zeros(N,Q);
        for i=1:Q
            r(:,i) = compute_radon_transform(f,theta(i),transposed);
        end
    else
        r = zeros(N,N);
        for i=1:Q
            r = r + compute_radon_transform(f(:,i),theta(i),transposed);
        end        
    end
    return;
end


%% basic operators %%

[Y,X] = meshgrid(1:N,1:N);
% indexing
V = @(theta)Y + (X-1)*tan(theta);
V = @(theta)mod(round(V(theta))-1,N)+1;
flat = @(x)x(:);
U = @(theta)reshape( sub2ind( [N N], X(:), flat(V(theta))), [N N]);
% Radon transform along one direction
R = @(f,theta)sum(f(U(theta)))';
% Adjoint operator
Rs = @(r,theta)assign_val(zeros(N), U(theta), repmat( r(:)', [N 1] ) ); 
reverse = @(x)x(end:-1:1);

% normalize the theta in [-pi/4,]
theta = mod(theta,pi);
if theta>3*pi/4
    theta = theta-pi;
end

if transposed<=0
    % forward
    if theta<=pi/4
        r = R(f, theta);
    else % if theta<=pi/2
        r = reverse(R(f', pi/2-theta));
    end
else
    % transpose
    if theta<=pi/4
        r = Rs(f, theta);
    else
        r = Rs(reverse(f) , pi/2-theta)';
    end
end