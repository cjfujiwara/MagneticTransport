%% Helmholtz

C1=[0 -.5 1];
C2=[0 .5 1];

[x,z]=meshgrid(linspace(-1.25,1.25,30),linspace(-1.25,1.25,30));
[Bx1,Bz1]=fieldCoil_2D(x,z,C1);
[Bx2,Bz2]=fieldCoil_2D(x,z,C2);

[~,B01]=fieldCoil_2D(0,0,C1);
[~,B02]=fieldCoil_2D(0,0,C2);

B0=B01+B02;

levels=linspace(0,1,50)*B0;

Bx=Bx1+Bx2;
Bz=Bz1+Bz2;

% 
zlin=linspace(-2,2,100);
[Bx1lin,Bz1lin]=fieldCoil_2D(zlin*0,zlin,C1);
[Bx2lin,Bz2lin]=fieldCoil_2D(zlin*0,zlin,C2);
Bxlin=Bx1lin+Bx2lin;
Bzlin=Bz1lin+Bz2lin;


figure(2)
clf

subplot(221);
quiver(x,z,Bx,Bz)


subplot(222);

plot(zlin,Bzlin);

subplot(223);
B=sqrt(Bx.^2+Bz.^2);

contour(x,z,B,100)
