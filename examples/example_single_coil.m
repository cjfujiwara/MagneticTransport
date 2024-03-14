
% 1 m radius coil at the origin
C=[0 0 1];

% Calculate the magnetic field in X-Z plane
[x,z]=meshgrid(linspace(-1.25,1.25,50),linspace(-.5,.5,50));
[Bx,Bz]=fieldCoil_2D(x,z,C);
B=sqrt(Bx.^2+Bz.^2);

% Calculate magnetic field along z axis
zlin=linspace(-2,2,100);
[Bxlin,Bzlin]=fieldCoil_2D(zlin*0,zlin,C);

% Calculate magnetic field at origin
[~,B0]=fieldCoil_2D(0,0,C);
levels=linspace(0,1,50)*B0;

% Make a plot
figure(1)
clf
subplot(221);
quiver(1e2*x,1e2*z,100*Bx,100*Bz)
xlabel('position (cm)')
ylabel('position (cm)')

subplot(222);
plot(zlin*1e2,Bzlin*1e2);
xlabel('axial position (cm)')
ylabel('magnetic field (G)')

subplot(223);
contour(x*1e2,z*1e2,B*1e2,100)
xlabel('position (cm)')
ylabel('position (cm)')
