
%% Anti Helmholtz

C1=[0 -.5 1];
C2=[0 .5 1];



% 
zlin=linspace(-2,2,100);
[Bx1lin,Bz1lin]=fieldCoil_2D(zlin*0,zlin,C1);
[Bx2lin,Bz2lin]=fieldCoil_2D(zlin*0,zlin,C2);
Bxlin=Bx1lin-Bx2lin;
Bzlin=Bz1lin-Bz2lin;


xlin=linspace(-2,2,100);
[Bx1lin2,Bz1lin2]=fieldCoil_2D(xlin,xlin*0,C1);
[Bx2lin2,Bz2lin2]=fieldCoil_2D(xlin,xlin*0,C2);
Bxlin2=Bx1lin2-Bx2lin2;
Bzlin2=Bz1lin2-Bz2lin2;

figure(3)
clf
set(gcf,'color','w');

subplot(221);
[x,z]=meshgrid(linspace(-1.25,1.25,16),linspace(-1.25,1.25,16));
[Bx1,Bz1]=fieldCoil_2D(x,z,C1);
[Bx2,Bz2]=fieldCoil_2D(x,z,C2);
Bx=Bx1-Bx2;
Bz=Bz1-Bz2;
quiver(x,z,Bx,Bz)
xlabel('x')
ylabel('z');
set(gca,'fontsize',14,'box','on','linewidth',1,'fontname','times');

subplot(222);
plot(zlin,Bzlin,'linewidth',1);
xlabel('z')
ylabel('Bz');
set(gca,'fontsize',14,'box','on','linewidth',1,'fontname','times');

subplot(224);
plot(xlin,Bxlin2,'linewidth',1);
xlabel('x')
ylabel('Bx');
set(gca,'fontsize',14,'box','on','linewidth',1,'fontname','times');


subplot(223);
[x,z]=meshgrid(linspace(-1.25,1.25,100),linspace(-1.25,1.25,100));
[Bx1,Bz1]=fieldCoil_2D(x,z,C1);
[Bx2,Bz2]=fieldCoil_2D(x,z,C2);
Bx=Bx1-Bx2;
Bz=Bz1-Bz2;
B=sqrt(Bx.^2+Bz.^2);
contour(x,z,B,500)
xlabel('x')
ylabel('z');
set(gca,'fontsize',14,'box','on','linewidth',1,'fontname','times');

