% example_horizontal_transport.m
%
% Horizontal magnetic transport is acheived by successively ramping the
% current through three pairs of anti-helmholtz coils.  As with all
% transport, we wish to smooth move the magnetic field zero, while keeping
% a fixed gradient in both the axial and transvese directions.


coils = makeHorizontalCoils;

ind0 = 2;

inds = [ind0 ind0+1 ind0+2];


C1 = coils(inds(1)).Coil;
C2 = coils(inds(2)).Coil;
C3 = coils(inds(3)).Coil;
x1=0.5*(coils(inds(1)).Position(1)+coils(inds(2)).Position(1));
x2=0.5*(coils(inds(2)).Position(1)+coils(inds(3)).Position(1));

i1 = 1;
i2 = 2;
i3 = 3;
s1 = coils(inds(1)).Name;
s2 = coils(inds(2)).Name;
s3 = coils(inds(3)).Name;

xvec = linspace(x1,x2,100);

%% Solving Tranport

imat = zeros(3,length(xvec));

% xvec=x1;
for jj = 1:length(xvec)
    z0 = 0;
    x0 = xvec(jj);
    
    % Distance to calculate taylor expansion
    dL  = 0.001;
    
    A = linspace(0,1,100);
    
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
    
    M1 = [Bx1 Bx2 Bx3;
        By1 By2 By3;
        Bz1 Bz2 Bz3];
    M1 = M1*1e4; % Convert Tesla to Gauss
    
    M2 = [Gx1 Gx2 Gx3;
        Gy1 Gy2 Gy3;
        Gz1 Gz2 Gz3];
    M2 = M2*100; % Convert Tesla/m to Gauss/cm

    % Target 100 G/cm in z direction
    Gz = 100;

    % Ratio of Gx/Gy
    % alpha=1/3; 

    alpha=1/1.705; 
    % alpha=1/1.0605;
    % alpha=1/1.45;

    % Calculate other gradients via Gauss's law
    Gy = -Gz/(1+alpha);
    Gx = alpha*Gy;

    b = [0; 0; 0; Gx; Gy; Gz];
    Amat = [M1;M2];
    isolve= mldivide(Amat,b);

    imat(:,jj)=isolve;
end

figure(6);
clf

co=get(gca,'colororder');
plot(xvec,imat(1,:),'color',co(i1,:)); hold on
plot(xvec,imat(2,:),'color',co(i2,:)); hold on
plot(xvec,imat(3,:),'color',co(i3,:)); hold on
xlabel('position (mm)')
xlim([x1 x2]);
ylabel('current (A)')
ylim([0 150]);
legend({s1,s2,s3})

