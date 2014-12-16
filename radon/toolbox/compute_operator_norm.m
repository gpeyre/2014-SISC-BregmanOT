function [L,e] = compute_operator_norm(A,n, niter)

% compute_operator_norm - compute operator norm
%
%   [L,e] = compute_operator_norm(A,n);
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<3
    niter = 30;
end

if length(n)==1
    u = randn(n,1); u = u/norm(u(:));
else
    u = n;
    u = u/norm(u(:));
end
e = [];
for i=1:30
    v = A(u);
    e(end+1) = sum(u(:).*v(:));
    u = v/norm(v(:));
end
L = e(end);