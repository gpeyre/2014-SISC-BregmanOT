epsilon=1.e-3;
N=100;
mass=0.7;
eta=1.e-4;
tol=1.e-4;
[x1,x2,y1,y2,xI,yI,xf,yf,mmP,m1,n1,totalTime,iteration]= partialOT(tol,epsilon,N,mass,eta);

str = sprintf('Total Time: %d',  totalTime);
disp(str);
str = sprintf('Iterations: %d',  iteration);
disp(str);

plot(x1,x2,'.k',y1,y2,'.k',xI,yI,'.r',xf,yf,'.g')

%write a file .png with the active and the inactive regions
%imm=mat2gray(mmP');
%imm=(imm-max(max(imm)))*(-1); 
%[I,m]=gray2ind(imm,16); 
%map=zeros(16,3);
%map(1,:)=[1 0 0];
%map(6,:)=[0 1 0];
%map(11,:)=[0 0 0];
%map(16,:)=[1 1 1]; 
%imwrite(I,map,'ActiveRegionsPartial.png')
