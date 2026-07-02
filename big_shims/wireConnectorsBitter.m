%% Introduction
% Want to estimate the fringing fields induced by input/output connectors
% for a feshbach configured coil.
%
%
% 2026/07/01 
% Distance of closet approach is 2.5", the magnet thickness is 1.2", and
% the OD of the magnet is 7.88"
%

OD = 7.88 * 0.0254;
R = OD/2;
t = 1.2 * 0.0254;
D = 2.5 * 0.0254;

% assume that the atoms are at x=y=z=0
ratom = [0 0 0];
%% Simple Wires at y=0 plane.
% There are various geometries to consider. But for now consider just two
% wires that stick out radially in the x-z plane.

% Top Magnet, Closer Wire
r1a = ratoms + [R 0 0] + [0 0 D/2];                % Start
r1b = ratoms + [R 0 0] + [0 0 D/2] + [10 0 0];     % End

% Top Magnet, Further Wire
r2a = ratoms + [R 0 0] + [0 0 D/2+t];              % Start
r2b = ratoms + [R 0 0] + [0 0 D/2+t] + [10 0 0];   % End

% Bottom Magnet, Closer Wire
r3a = ratom + [R 0 0] + [0 0 -D/2];               % Start
r3b = ratom + [R 0 0] + [0 0 -D/2] + [10 0 0];    % End

% Bottom 
r4a = ratom + [R 0 0] + [0 0 -D/2-t];             % Start
r4b = ratom + [R 0 0] + [0 0 -D/2-t] + [10 0 0];  % End

%%
A = 1;

x = linspace(-1,9,1e3)*ft2m;
y = linspace(-1,6,1e3)*ft2m;
z = 16.5/12*ft2m;
[X,Y,Z] = meshgrid(x,y,z);
[Bx,By,Bz] = fieldLine(X,Y,Z,Lm,Rm,A);
Bx= Bx*1e3;
By= By*1e3;
Bz= Bz*1e3;

f=figure(10);
clf
f.Color='w';

subplot(211);
imagesc(x*12/ft2m,y*12/ft2m,Bz);
axis equal tight
cc=colorbar;
cc.Label.String = 'magnetic field (mG)';
caxis([0 10]);
xlabel('x position (in.)');
ylabel('y position (in.)');
title('Magnetic Field B_z @ 1A (mG); z=16.5 in.');
set(gca,'yDir','normal');

str = [num2str(l*12) ' in. x ' num2str(w*12) ' in. coils,' num2str(dZ*12) ' in. separation'];
tt = text(5,5,str,'units','pixels','fontsize',14,'color','w',...
    'horizontalalignment','left','verticalalignment','bottom');
hold on
plot(atom(1),atom(2),'ro','markerfacecolor','r','markersize',10);
plot(sensor(1),sensor(2),'bo','markerfacecolor','b','markersize',10);

subplot(234);
x = linspace(0,65,1e3)/12*ft2m;
y = w/2*ft2m;
z = 16.5/12*ft2m;
[Bx,By,Bz] = fieldLine(x,y,z,Lm,Rm,A);
Bx= Bx*1e3;
By= By*1e3;
Bz= Bz*1e3;
plot(12*x/ft2m,Bz)
xlabel('x position (in)');
ylabel('Bz (G)');
title(['B vs x at y=' num2str(y*12/ft2m) ' in, z=' num2str(z*12/ft2m) ' in']);
xlim([0 max(x*12/ft2m)]);
ylim([0 10]);

hold on
plot([1 1]*atom(1),get(gca,'YLim'),'r-');
plot([1 1]*sensor(1),get(gca,'YLim'),'b-');

subplot(235);
% x = linspace(0,65,1e3)/12*ft2m;
x = l/2*ft2m;
% y = w/2*ft2m;
y = linspace(0,w,1e3)*ft2m;
z = 16.5/12*ft2m;
[Bx,By,Bz] = fieldLine(x,y,z,Lm,Rm,A);
Bx= Bx*1e3;
By= By*1e3;
Bz= Bz*1e3;
plot(12*y/ft2m,Bz)
xlabel('x position (in)');
ylabel('Bz (G)');
title(['B vs x at x=' num2str(x*12/ft2m) ' in']);
xlim([0 max(y*12/ft2m)]);
ylim([0 10]);

subplot(236);
% x = linspace(0,65,1e3)/12*ft2m;
x = l/2*ft2m;
 y = w/2*ft2m;
% y = linspace(0,w,1e3)*ft2m;
% z = 16.5/12*ft2m;
z = linspace(10,16.5,1e3)/12*ft2m;

[Bx,By,Bz] = fieldLine(x,y,z,Lm,Rm,A);
Bx= Bx*1e3;
By= By*1e3;
Bz= Bz*1e3;
plot(12*z/ft2m,Bz)
xlabel('x position (in)');
ylabel('Bz (G)');
title(['B vs x at x=' num2str(x*12/ft2m) ' in']);
xlim([min(z) max(z)]*12/ft2m);
ylim([0 10]);

%%
x = 0;
y = linspace(-1,1,1e3);
z = 0;
A = 1;

[Bx,By,Bz] = fieldLine(x,y,z,L,R,A);

figure(3)
clf
% plot(y,Bx);
% hold on
% plot(y,By);
plot(y,Bz);
ylim([-.01 .01]);

%%

x = linspace(-201,201,1e3+1);
y = linspace(-10,10,100);
z = 0;
[X,Y,Z]=meshgrid(x,y,z);

[Bx,By,Bz] = fieldLine(X,Y,Z,L,R,A);

figure(2)
clf

subplot(121);
imagesc(x,y,Bz)
xlabel('x');
ylabel('y');

subplot(122);

plot(y,Bz(:,501))
hold on
plot(y,Bz(:,230))

% xlabel('x');
% ylabel('y');

%%

x = 0;
y = linspace(-1,1,100);
z = linspace(-1,1,100);
[X,Y,Z]=meshgrid(x,y,z);

[Bx,By,Bz] = fieldLine(X,Y,Z,L,R,A);

figure(3)
clf
quiver(Z,Y,Bz,By)
xlabel('z');
ylabel('y');
