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

Iamp = [1 1];

I1 = Iamp(1);
I2 = Iamp(2);

C1 = [0 1*R R 1];
C2 = [0 2*R R 1];
C3 = [0 3*R R 1];
C4 = [0 4*R R 1];

%% Solving Tranport
zvec = linspace(R,2*R,100);

imat = zeros(4,length(zvec));
for jj = 1:length(zvec)
    z0 = zvec(jj);
    x0 = 0;
    
    % Distance to calculate taylor expansion
    dL  = 0.001;
    
    A = linspace(0,1,100);
    
    [Bx1,By1,Bz1]=fieldCoil_3D(x0,0,z0,C1);
    [Bx2,By2,Bz2]=fieldCoil_3D(x0,0,z0,C2);
    [Bx3,By3,Bz3]=fieldCoil_3D(x0,0,z0,C3);
    [Bx4,By4,Bz4]=fieldCoil_3D(x0,0,z0,C4);

    % X Grad
    [Bx1p,~,~]=fieldCoil_3D(x0+dL,0,z0,C1);
    [Bx2p,~,~]=fieldCoil_3D(x0+dL,0,z0,C2);
    [Bx3p,~,~]=fieldCoil_3D(x0+dL,0,z0,C3);
    [Bx4p,~,~]=fieldCoil_3D(x0+dL,0,z0,C4);

    [Bx1n,~,~]=fieldCoil_3D(x0-dL,0,z0,C1);
    [Bx2n,~,~]=fieldCoil_3D(x0-dL,0,z0,C2);
    [Bx3n,~,~]=fieldCoil_3D(x0-dL,0,z0,C3);
    [Bx4n,~,~]=fieldCoil_3D(x0-dL,0,z0,C4);

    Gx1 = (Bx1p-Bx1n)/(2*dL);
    Gx2 = (Bx2p-Bx2n)/(2*dL);
    Gx3 = (Bx3p-Bx3n)/(2*dL);
    Gx4 = (Bx4p-Bx4n)/(2*dL);

    
    % Y Grad
    [~,By1p,~]=fieldCoil_3D(x0,dL,z0,C1);
    [~,By2p,~]=fieldCoil_3D(x0,dL,z0,C2);
    [~,By3p,~]=fieldCoil_3D(x0,dL,z0,C3);
    [~,By4p,~]=fieldCoil_3D(x0,dL,z0,C4);

    [~,By1n,~]=fieldCoil_3D(x0,-dL,z0,C1);
    [~,By2n,~]=fieldCoil_3D(x0,-dL,z0,C2);
    [~,By3n,~]=fieldCoil_3D(x0,-dL,z0,C3);
    [~,By4n,~]=fieldCoil_3D(x0,-dL,z0,C4);

    Gy1 = (By1p-By1n)/(2*dL);
    Gy2 = (By2p-By2n)/(2*dL);
    Gy3 = (By3p-By3n)/(2*dL);
    Gy4 = (By4p-By4n)/(2*dL);

    % Z Grad
    [~,~,Bz1p]=fieldCoil_3D(x0,0,z0+dL,C1);
    [~,~,Bz2p]=fieldCoil_3D(x0,0,z0+dL,C2);
    [~,~,Bz3p]=fieldCoil_3D(x0,0,z0+dL,C3);
    [~,~,Bz4p]=fieldCoil_3D(x0,0,z0+dL,C4);
    [~,~,Bz1n]=fieldCoil_3D(x0,0,z0-dL,C1);
    [~,~,Bz2n]=fieldCoil_3D(x0,0,z0-dL,C2);
    [~,~,Bz3n]=fieldCoil_3D(x0,0,z0-dL,C3);
    [~,~,Bz4n]=fieldCoil_3D(x0,0,z0-dL,C4);

    Gz1 = (Bz1p-Bz1n)/(2*dL);
    Gz2 = (Bz2p-Bz2n)/(2*dL);
    Gz3 = (Bz3p-Bz3n)/(2*dL);
    Gz4 = (Bz4p-Bz4n)/(2*dL);
    
    M1 = [Bx1 Bx2 Bx3 Bx4;
        By1 By2 By3 By4;
        Bz1 Bz2 Bz3 Bz4];
    
    M2 = [Gx1 Gx2 Gx3 Gx4;
        Gy1 Gy2 Gy3 Gy4;
        Gz1 Gz2 Gz3 Gz4];
    
    b = [0; 0; 0; 1; -.5; -.5];
    Amat = [M1;M2];
    isolve= mldivide(Amat,b);

    imat(:,jj)=isolve;
end

figure(5);
clf
plot(zvec,imat(1,:)); hold on
plot(zvec,imat(2,:)); hold on
plot(zvec,imat(3,:)); hold on
plot(zvec,imat(4,:)); hold on
