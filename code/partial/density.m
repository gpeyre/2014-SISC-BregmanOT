function [mu,nu,areaX,areaY]=density(dx,dy)
%density
% dx: displacement for the first ellipse
% dy: displacement for the second ellipse


%ellipses
Mx = [0.2 -0.1; 0 0.1];
Ms=[0.2 0;0 0.2 ];     %use Ms to create an annulus
My =  [0.25 0; 0 0.25];
Ms1= [0.15 0; 0 0.15];
Qx = inv(Mx);
Qy = inv(My);
Qs=inv(Ms);
Qs1=inv(Ms1);
 a = Qx(1,1)^2+Qx(2,1)^2; d = Qx(1,1)*Qx(1,2)+Qx(2,1)*Qx(2,2); c = Qx(1,2)^2+Qx(2,2)^2;
 a1 = Qs(1,1)^2+Qs(2,1)^2; d1 = Qs(1,1)*Qs(1,2)+Qs(2,1)*Qs(2,2); c1 = Qs(1,2)^2+Qs(2,2)^2;
 areaX = pi/sqrt(a*c-d^2);%-pi/sqrt(a1*c1-d1^2);
 mu = @(x,y) (((Qx(1,1)*(x-0.25)+Qx(1,2)*(y-0.25)).^2+(Qx(2,1)*(x-0.25)+Qx(2,2)*(y-0.25)).^2)<=1);%-(((Qs(1,1)*x+Qs(1,2)*y).^2+(Qs(2,1)*x+Qs(2,2)*y).^2)<=1);
% ms = @(x,y) (((Qs(1,1)*(x-0.2)+Qs(1,2)*(y-0.2)).^2+(Qs(2,1)*(x-0.2)+Qs(2,2)*(y-0.2)).^2)<=1);

 a = Qy(1,1)^2+Qy(2,1)^2; d = Qy(1,1)*Qy(1,2)+Qy(2,1)*Qy(2,2); c = Qy(1,2)^2+Qy(2,2)^2;
 a1 = Qs1(1,1)^2+Qs1(2,1)^2; d1 = Qs1(1,1)*Qs1(1,2)+Qs1(2,1)*Qs1(2,2); c1 = Qs1(1,2)^2+Qs1(2,2)^2;
 areaY = pi/sqrt(a*c-d^2);%-pi/sqrt(a1*c1-d1^2);

 nu = @(x,y) (((Qy(1,1)*(x-dx)+Qy(1,2)*(y-dy)).^2+(Qy(2,1)*(x-dx)+Qy(2,2)*(y-dy)).^2)<=1);%-(((Qs1(1,1)*(x-dx)+Qs1(1,2)*(y-dy)).^2+(Qs1(2,1)*(x-dx)+Qs1(2,2)*(y-dy)).^2)<=1);


 %for the diamond or the square
 %mu=@(x,y) (x>=0 & x<=1/6).*(y>=0 & y<=0.5)+(y>=0 & y<=1/6).*(x>=0 & x<=1/3);
 %areaX=1/9;
 %nu=@(x,y) (x>=0.5 & x<=2/3).*(y>=0.5 & y<=1.0)+(x>=5/6 & x<=1).*(y>=0.5 & y<=1)+(x>2/3 & x<5/6).*(y>=-x+3/2 & y<=-x+5/3);
