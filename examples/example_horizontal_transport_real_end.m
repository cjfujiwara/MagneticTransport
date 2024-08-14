% example_horizontal_transport.m
%
% Horizontal magnetic transport is acheived by successively ramping the
% current through three pairs of anti-helmholtz coils.  As with all
% transport, we wish to smooth move the magnetic field zero, while keeping
% a fixed gradient in both the axial and transvese directions.


[~,coils] = makeHorizontalCoils;


C1 = coils(end-1).Coil;
C2 = coils(end).Coil;


x1=0.5*(coils(inds(end-1)).Position(1)+coils(inds(end)).Position(1));
x2=coils(inds(end)).Position(1);

i1 = 1;
i2 = 2;
s1 = coils(inds(end-1)).Name;
s2 = coils(inds(end)).Name;

xvec = linspace(x1,x2,100);

%% Solving Tranport

imat = zeros(2,length(xvec));
Gx = zeros(1,length(xvec));
Gy = zeros(1,length(xvec));
Gz = zeros(1,length(xvec));

% xvec=x1;
for jj = 1:length(xvec)
    z0 = 0;
    x0 = xvec(jj);
    
    % Distance to calculate taylor expansion
    dL  = 0.001;
    
    A = linspace(0,1,100);
    
    [Bx1,By1,Bz1]=fieldCoil_3D(x0,0,z0,C1);
    [Bx2,By2,Bz2]=fieldCoil_3D(x0,0,z0,C2);
    
    % X Grad
    [Bx1p,~,~]=fieldCoil_3D(x0+dL,0,z0,C1);
    [Bx2p,~,~]=fieldCoil_3D(x0+dL,0,z0,C2);
    [Bx1n,~,~]=fieldCoil_3D(x0-dL,0,z0,C1);
    [Bx2n,~,~]=fieldCoil_3D(x0-dL,0,z0,C2);
    Gx1 = 100*(Bx1p-Bx1n)/(2*dL);
    Gx2 = 100*(Bx2p-Bx2n)/(2*dL);
    
    % Y Grad
    [~,By1p,~]=fieldCoil_3D(x0,dL,z0,C1);
    [~,By2p,~]=fieldCoil_3D(x0,dL,z0,C2);
    [~,By1n,~]=fieldCoil_3D(x0,-dL,z0,C1);
    [~,By2n,~]=fieldCoil_3D(x0,-dL,z0,C2);
    Gy1 = 100*(By1p-By1n)/(2*dL);
    Gy2 = 100*(By2p-By2n)/(2*dL);
    
    % Z Grad
    [~,~,Bz1p]=fieldCoil_3D(x0,0,z0+dL,C1);
    [~,~,Bz2p]=fieldCoil_3D(x0,0,z0+dL,C2);
    [~,~,Bz1n]=fieldCoil_3D(x0,0,z0-dL,C1);
    [~,~,Bz2n]=fieldCoil_3D(x0,0,z0-dL,C2);
    Gz1 = 100*(Bz1p-Bz1n)/(2*dL);
    Gz2 = 100*(Bz2p-Bz2n)/(2*dL);
    
    M1 = [Bx1 Bx2;
        By1 By2;
        Bz1 Bz2];
    M1 = M1*1e4; % Convert Tesla to Gauss
    
    M2 = [Gz1 Gz2];
    % M2 = M2*100; % Convert Tesla/m to Gauss/cm

    % Target 100 G/cm in z direction
    Gz0 = 100;

    % Ratio of Gx/Gy
    % alpha=1/1.705; 
    % alpha=1/1.0605;
    % alpha=1/1.45;

    % Calculate other gradients via Gauss's law
    % Gy = -Gz/(1+alpha);
    % Gx = alpha*Gy;

    b = [0; 0; 0; Gz0];
    Amat = [M1;M2];
    isolve= mldivide(Amat,b);

    Gx(jj) = isolve(1)*Gx1+isolve(2)*Gx2;
    Gy(jj) = isolve(1)*Gy1+isolve(2)*Gy2;
    Gz(jj) = isolve(1)*Gz1+isolve(2)*Gz2;

    imat(:,jj)=isolve;
end

figure(7);
clf

subplot(211);
co=get(gca,'colororder');
plot(xvec,imat(1,:),'color',co(i1,:)); hold on
plot(xvec,imat(2,:),'color',co(i2,:)); hold on
% plot(xvec,imat(3,:),'color',co(i3,:)); hold on
xlabel('position (mm)')
xlim([x1 x2]);
ylabel('current (A)')
ylim([0 150]);
legend({s1,s2,s3})

subplot(212);
co=get(gca,'colororder');
plot(xvec,Gy./Gx,'color',co(i1,:)); hold on
% plot(xvec,Gy,'color',co(i2,:)); hold on
% plot(xvec,Gz,'color',co(i3,:)); hold on


