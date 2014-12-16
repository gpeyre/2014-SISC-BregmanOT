rand('seed',0);
close all

N=30;
NN=60;

xyes=1./(1+exp(5*(rand(1,N)-.5)));
yyes=1./(1+exp(5*(rand(1,N)-.5)));

sizeScale=.2*rand(N,1)+.8;

imgNNNN=zeros(NN,NN,N);
factor=2.3;
figure(1);

myColorMap = hot; % Make a copy of jet.
% Assign white (all 1's) to black (the first row in myColorMap).
%myColorMap(1, :) = [1 1 1];
myColorMap=flipud(myColorMap);

for i=1:N,
    a1=factor*(8*rand+5);
    b1=factor*(8*rand+5);
    angle=360*rand;
    p1 = calculateEllipse(0, 0, a1, b1, angle ,500);
    %p1 = [p1;calculateEllipse(0, 0, .96*a1, .96*b1, angle ,500)];
    
    a3=factor*(5*rand+4);
    b3=factor*(5*rand+4);
    angle=360*rand;
    x=rand();
    y=rand();
    scalx=(rand()-.5)*.8;
    scaly=(rand()-.5)*.8;
    p3= calculateEllipse(scalx*a1 , scaly*b1, a3, b3, angle ,500);
    %p3= [p3;calculateEllipse(scalx*a1 , scaly*b1, .93*a3, .93*b3, angle ,500)];
    
    
    a=factor*(2*rand+2);
    b=factor*(2*rand+2);
    angle=360*rand;
    x=rand();
    y=rand();
    scalx=(rand()-.5)*.8;
    scaly=(rand()-.5)*.8;
    p2= calculateEllipse(scalx*a3 , scaly*b3, a, b, angle ,500);
    %p2= [p2;calculateEllipse(scalx*a3 , scaly*b3, .93*a, .93*b, angle ,500)]; % add more points to have smoother looking ellipse
    
    
    xoffset=min([p1(:,1);p2(:,1);p3(:,1)]);
    yoffset=min([p1(:,2);p2(:,2);p3(:,2)]);
    p1(:,1)=p1(:,1)-xoffset+1;
    p2(:,1)=p2(:,1)-xoffset+1;
    p3(:,1)=p3(:,1)-xoffset+1;
    
    p1(:,2)=p1(:,2)-yoffset+1;
    p2(:,2)=p2(:,2)-yoffset+1;
    p3(:,2)=p3(:,2)-yoffset+1;
    
    
    
    %     plot(p1(:,1), p1(:,2), '.-'), axis equal;
    %     hold on;
    %     plot(p2(:,1), p2(:,2), '.-'), axis equal;
    %     pause();
    p1=round(p1); p2=round(p2); p3=round(p3);
    
    
    xmax=max([p1(:,1);p2(:,1);p3(:,1)]);
    ymax=max([p1(:,2);p2(:,2);p3(:,2)]);
    
    IMG=zeros(xmax,ymax);
    for j=1:size(p1,1),
        IMG(p1(j,1),p1(j,2))=1;
        IMG(p2(j,1),p2(j,2))=1;
        IMG(p3(j,1),p3(j,2))=1;
    end
    IMG=IMG/sum(sum(IMG));
    
    %     xyes=rand()%>.5;
    %     yyes=rand()%>.5;
    Im=imresize(IMG,sizeScale(i));
    %Im=IMG;
    Im(Im<1e-15)=0;
    [nIm,mIm]=size(Im);
    
    startn=round((NN-nIm)*xyes(i));
    startm=round((NN-mIm)*yyes(i));
    imgNNNN(1+startn:nIm+startn,1+startm:mIm+startm,i)=Im;
    imgNNNN(:,:,i)=imgNNNN(:,:,i)/sum(sum(imgNNNN(:,:,i)));
    subaxis(6,10,floor((i-1)/6)+1,mod(i-1,6)+1,'Spacing', 0.005, 'Padding', 0, 'Margin', 0.005);
    
    imagesc(imgNNNN(:,:,i));
    colormap(myColorMap);    
    box on;
    set(gca,'xtick',[],'ytick',[])
    xlabel('');
    ylabel('');
    
    pause(.015);
end

R=reshape(imgNNNN,NN^2,N);

%NMAX=2000;
%gpuDevice( 3 )
%reset(gpuDevice( 1 ))

M=computeDistanceMatrixGrid(NN,NN).^2;
M=M/median(M(:));

toleranceDifference=1e-4;
inLoopTolerance=.1;

lambda=1000;

% try gpuDevice(1),
%     USINGGPU=true
% catch ME
%     USINGGPU=false
% end



StartingTrick=1;

t0=1;
%c1=smoothConjugateDiscreteWassersteinBarycenter(R,M,100000,t0,lambda,USINGGPU,toleranceDifference,StartingTrick);
cb=bregmanWassersteinBarycenter(R,M,100000,lambda,USINGGPU,toleranceDifference*.0001);
%c0=smoothDualPrimalWassersteinBarycenter(R,M,100000,lambda,USINGGPU,toleranceDifference);

%c2=discreteWassersteinBarycenter(R,M,lambda,50,t0,USINGGPU,toleranceDifference,inLoopTolerance);
%if USINGGPU,
%    c2=gather(c2);
%end
%c3=discreteWassersteinBarycenter(R,M,100,.1);

%c=gather(acceleratedDiscreteWassersteinBarycenter(R,M,lambda,50,t0,USINGGPU));
