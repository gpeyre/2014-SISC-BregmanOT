%%
% Test for the construction from a discrete Radon transform
% using Wasserstein fidelity.

addpath('toolbox/');

% size of the image
N = 128/2; 
N = 100;
N = 60;

%% useful helpers %%
% remove ticks
remove_ticks = @()set(gca, 'XTick', [], 'YTick', []);
% set aspect ratio
set_ar = @()set(gca, 'PlotBoxAspectRatio', [1 1/3 1]);


% helpers
dotp = @(x,y)sum(x(:).*y(:));
flat = @(x)x(:);
mynorm = @(x)norm(x(:));
R  = @(f,theta)compute_radon_transform(f,theta,0);
Rs = @(r,theta)compute_radon_transform(r,theta,1);
% set figure title
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

% projection on pi*1=p
ProjR = @(pi,p)pi .* repmat( p./max(sum(pi,2), 1e-10), [1 size(pi,2)] );
% projection on pi^T*1=q
ProjC = @(pi,q)pi .* repmat( q'./max(sum(pi,1), 1e-10), [size(pi,1) 1] );

name = 'disk';
name = '2disks';
name = 'phantom';

rep = ['results/' name '/'];
if not(exist(rep))
    mkdir(rep);
end


switch name
    case 'disk'
        options.radius = .3;
        options.center = [.4 .45];
        f = rescale( load_image(name, N,options) );
        % template for the matching
        delta = round([.2 .2]*N);
        f0 = circshift(f, delta );
    case 'phantom'
        f = rescale( load_image(name, N) );
        % template for the matching
        options.radius = .4;
        options.center = [.5 .5];
        f0 = rescale( load_image('disk', N,options) );
end
f = f/sum(f(:));
f0 = f0/sum(f0(:));

if 0
    % do some padding
    P = ceil(N/16)*2; P=0;
    f0 = load_image(name, N-P,options);
    f = zeros(N);
    f(P/2+1:end-P/2,P/2+1:end-P/2) = f0;
end
        
%%
% Display of a dense Radon transform

Q = 128;
Theta = linspace(-pi/4,3*pi/4,Q+1); Theta(1) = [];

clf;
imageplot(R(f,Theta));
drawnow;

% check for adjoint correctness
g = randn(N); r = randn(N,length(Theta));
e = dotp( R(g,Theta),r ) - dotp(g,Rs(r,Theta) );
if abs(e)>1e-9
    warning('Problem with adjointness');
end


%%
% The goal is to solve, given some input {r_k}_k
%   min_f lambda*W(f,f0) + sum_k W(R(f,theta_k),r_k)
% which we rescast as
%   min_{Pi, pi_theta}  lambda*KL(Pi,bar Pi) + sum_k KL(pi_k,bar pi_k)
%       Pi'*1 = f0,     pi_k'*1 = r_k
%       pi_k*1 = R(f,theta_k)   and   Pi*1=f
%       ==>  pi_k*1=R(Pi*1,theta_k)

Q = 12;
Theta = linspace(-pi/4,3*pi/4,Q+1); Theta(1) = [];

% observations
rk = R(f,Theta);

%%
% Do a  L^2 reconstruction, not exactly the pseudo inverse actually

[L,e] = compute_operator_norm(@(f)Rs(R(f,Theta),Theta),randn(N,N), 20);
niter = 50;
% min |R(f)-rk|^2
f1 = zeros(N);
tau = 1/L;
for i=1:niter
    f1 = f1 - tau*Rs(R(f1,Theta)-rk,Theta);
end

%%
% Wasserstein reconstruction.

% fidelity weight
if not(exist('lambda'))
    lambda = 1;
end

% regularization weight
gamma = (4/N)^2;
gamma = (2/N)^2;

% Wasserstein metric exponent
Wexp = 2;

% 2D ground metric
t = linspace(0,1,N);
[X1,Y1,X2,Y2] = ndgrid(t,t,t,t);
M2 = sqrt( (X1-X2).^2 + (Y1-Y2).^2 ).^Wexp;
M2 = reshape(M2, [N^2 N^2]);
% 1D ground metrix, modulo N (on the circle)
[Y,X] = meshgrid(t,t);
M1 = abs(X-X');
M1(M1>1/2) = 1-M1(M1>1/2);
M1 = M1.^Wexp;
% initialization bar Pi and bar pi_k
Pi  = exp(-M2/gamma);
pik = repmat( exp(-M1/gamma), [1 1 Q] );

ErrPi = [];
Errpi = [];
Err = [];
niter = 5000;
for i=1:niter
    progressbar(i,niter);
    % independent projection on each marginals of columns
    % Pi'*1 = f0, pi_k'*1 = r_k
    ErrPi(i) = mynorm( sum(Pi)'-f0(:) );
    Pi = ProjC(Pi,f0(:));
    for k=1:Q
        Errpi(i,k) = mynorm( sum(pik(:,:,k))'-rk(:,k) );
        pik(:,:,k) = ProjC( pik(:,:,k), rk(:,k) ) ;
    end    
    % (lambda,1)-KL projection on joint constraint    
    % pi_k*1=R(Pi*1,theta_k)
    for k=1:Q
        % R(Pi*1,theta_k)
        a = R( reshape(sum(Pi'),[N N]), Theta(k) );
        % pi_k*1
        b = sum( pik(:,:,k)' )';
        % record error for convergence study
        Err(i,k) = mynorm(a-b);
        % weighted geometric mean
        c = exp( ( lambda*log(a) + log(b) )/(1+lambda) );
        % impose on pi_k
        pik(:,:,k) = ProjR( pik(:,:,k), c );
        % impose on Pi % TOCORRECT -- BUG HERE %
        u = Rs(c./max(a,1e-10),Theta(k));
        Pi = Pi .* repmat( u(:), [1 N^2] );
    end
end

%% 
% Displays.

% error display
figure(1); setfigname('Error decay');
clf;
subplot(3,1,1);
plot(log10(ErrPi)); axis tight;
subplot(3,1,2);
plot(log10(Errpi)); axis tight;
subplot(3,1,3);
plot(log10(Err)); axis tight;

% display result
figure(2); setfigname('Measure');
g = reshape(sum(Pi'),[N N]);
r = R(g, Theta);
imageplot(g);

% display marginals
k = round(Q/4);
figure(3); setfigname('Marginals');
clf; hold on;
plot(r(:,1:k:end), '-');
plot(rk(:,1:k:end), '--');
axis tight; box on;

% export marginals
for i=1:Q
    clf;hold on;
    plot(r(:,i), 'k-');
    plot(rk(:,i), 'k--');
    axis tight; box on; 
    remove_ticks();
    set_ar();
    saveas(gcf, [rep 'marginal-' num2str(i) '.eps'], 'eps');
end


%%
% Save results

clamp = @(x)max(0,min(1,x));

imwrite(rescale(f), [rep 'original.png']);
imwrite(rescale(f0), [rep 'template.png']);
imwrite(rescale(f1), [rep 'recovered-Q' num2str(Q) '-L2.png']);
imwrite(rescale(g), [rep 'recovered-Q' num2str(Q) '-W' num2str(Wexp) '-lambda' num2str(round(lambda*10)) '.png']);
