% example_horizontal_transport.m
%
% Horizontal magnetic transport is acheived by successively ramping the
% current through three pairs of anti-helmholtz coils.  As with all
% transport, we wish to smooth move the magnetic field zero, while keeping
% a fixed gradient in both the axial and transvese directions.

%% Specify coils parameters
% Define three pairs of 
% Coil radius of 10 cm
R=0.1;

Iamp = [1 1 1];

I1 = Iamp(1);
I2 = Iamp(2);
I3 = Iamp(3);

dZ =  R*0;
cc = 1;
% Pair 1 (x0,z0,R)
C1 = [-R -R/2 R*cc I1;
    -R R/2 R*cc -I1];
% Pair 2 (x0,z0,R)
C2=[0 -(R/2-dZ) R*cc I2;
    0 (R/2-dZ) R*cc -I2];
% Pair 3 (x0,z0,R)
C3=[R -R/2 R*cc I3;
    R R/2 R*cc -I3];


%% Define the position mesh

x = linspace(-4*R,4*R,50);
z = linspace(-R,R,50);

[xx,zz]=meshgrid(x,z);

%% Calculate the magnetic field of each coil

[Bx1,By1,Bz1]=fieldCoil_3D(xx,0,zz,C1);
[Bx2,By2,Bz2]=fieldCoil_3D(xx,0,zz,C2);
[Bx3,By3,Bz3]=fieldCoil_3D(xx,0,zz,C3);

Bx = Bx1+Bx2+Bx3;
By = By1+By2+By3;
Bz = Bz1+Bz2+Bz3;

figure(1)
clf
quiver(1e2*x,1e2*z,100*Bx,100*Bz)
xlabel('position (cm)')
ylabel('position (cm)')

%% Solving Tranport
xvec = linspace(-R,R,100);
% % xvec = -R/2;
% xvec=0;
% xvec=0;
imat = zeros(3,length(xvec));
for jj = 1:length(xvec)
    z0 = 0;
    x0 = xvec(jj);
    
    % Distance to calculate taylor expansion
    dL  = 0.00001;
        
    [Bx1,By1,Bz1]=fieldCoil_3D(x0,0,z0,C1);
    [Bx2,By2,Bz2]=fieldCoil_3D(x0,0,z0,C2);
    [Bx3,By3,Bz3]=fieldCoil_3D(x0,0,z0,C3);
    
    % X Grad
    [Bx1p,~,~]=fieldCoil_3D(x0+dL,0,z0,C1);
    [Bx2p,~,~]=fieldCoil_3D(x0+dL,0,z0,C2);
    [Bx3p,~,~]=fieldCoil_3D(x0+dL,0,z0,C3);
    [Bx1n,~,~]=fieldCoil_3D(x0-dL,0,z0,C1);
    [Bx2n,~,~]=fieldCoil_3D(x0-dL,0,z0,C2);
    [Bx3n,~,~]=fieldCoil_3D(x0-dL,0,z0,C3);
    Gx1 = (Bx1p-Bx1n)/(2*dL);
    Gx2 = (Bx2p-Bx2n)/(2*dL);
    Gx3 = (Bx3p-Bx3n)/(2*dL);
    
    % Y Grad
    [~,By1p,~]=fieldCoil_3D(x0,dL,z0,C1);
    [~,By2p,~]=fieldCoil_3D(x0,dL,z0,C2);
    [~,By3p,~]=fieldCoil_3D(x0,dL,z0,C3);
    [~,By1n,~]=fieldCoil_3D(x0,-dL,z0,C1);
    [~,By2n,~]=fieldCoil_3D(x0,-dL,z0,C2);
    [~,By3n,~]=fieldCoil_3D(x0,-dL,z0,C3);
    Gy1 = (By1p-By1n)/(2*dL);
    Gy2 = (By2p-By2n)/(2*dL);
    Gy3 = (By3p-By3n)/(2*dL);
    
    % Z Grad
    [~,~,Bz1p]=fieldCoil_3D(x0,0,z0+dL,C1);
    [~,~,Bz2p]=fieldCoil_3D(x0,0,z0+dL,C2);
    [~,~,Bz3p]=fieldCoil_3D(x0,0,z0+dL,C3);
    [~,~,Bz1n]=fieldCoil_3D(x0,0,z0-dL,C1);
    [~,~,Bz2n]=fieldCoil_3D(x0,0,z0-dL,C2);
    [~,~,Bz3n]=fieldCoil_3D(x0,0,z0-dL,C3);
    Gz1 = (Bz1p-Bz1n)/(2*dL);
    Gz2 = (Bz2p-Bz2n)/(2*dL);
    Gz3 = (Bz3p-Bz3n)/(2*dL);

    G1 = [Gx1;Gy1;Gz1];
    G2 = [Gx2;Gy2;Gz2];
    G3 = [Gx3;Gy3;Gz3];

    B1 = [Bx1;By1;Bz1];
    B2 = [Bx2;By2;Bz2];
    B3 = [Bx1;By3;Bz3];

    M1 = [Bx1 Bx2 Bx3;
        By1 By2 By3;
        Bz1 Bz2 Bz3];
    M1 = M1*1e4;
    
    M2 = [Gx1 Gx2 Gx3;
        Gy1 Gy2 Gy3;
        Gz1 Gz2 Gz3];
    M2 = M2*1e2;


    % Target 100 G/cm in z direction
    Gz = 100;

    % Ratio of Gx/Gy
    alpha=1/1.665; 

    alpha=1/1.2; 

    % alpha=5;

    % Calculate other gradients via Gauss's law
    Gy = -Gz/(1+alpha);
    Gx = alpha*Gy;

    b = [0; 0; 0; Gx; Gy; Gz];
    Amat = [M1;M2];
    isolve= mldivide(Amat,b);

    imat(:,jj)=isolve;
end

figure(5);
clf
plot(xvec,imat(1,:)); hold on
plot(xvec,-imat(2,:)); hold on
plot(xvec,imat(3,:)); hold on

%% 

A = linspace(0,1,1000);

B = zeros(6,length(A));
for kk=1:length(A)
    b=[M1;M2]*[-A(kk); 1; -A(kk)];
    B(:,kk) = b;
end

