% Add to Path

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));

%% Construct the coils

% Solve vertical transport in the "easy" regime
coils = makeHorizontalCoils;
warning off;

options = optimoptions(@fmincon,'MaxIterations',5e4,'MaxFunctionEvaluations',1e5);

%% Establish Symmetry Points

x0 = 0;0;        % MOT location
xa = 0.04;0.0455;.044;       % [x0,xa] : [P,MOT,3]; to xa
xb = 0.0681;     % [xa,xb] : [MOT,3,4]; to 3/4 center
xc = 0.1006;     % [xb,xc] : [3,4,5]; to 4/5 center
xd = 0.1321;     % [xc,xd] : [4,5,6]; to 5/6 center
xe = 0.1636;     % [xd,xe] : [5,6,7]; to 6/7 center
xf = 0.1951;     % [xe,xf] : [6,7,8]; to 7/8 center
xg = 0.2266;     % [xf,xg] : [7,8,9]; to 8/9 center
xh = 0.2581;     % [xg,xh] : [8,9,10]; to 9/10 center
xi = 0.291;0.2959;     % [xh,xi] : [9,10,11]; to 10/11 cetner
xj = 0.3116;     % [xi,xj] : [10,11,11_ex]; to 11/11_ex cetner
xk = 0.339;      % [xj,xk] : [11,11_ex,12]; to 11_ex/12
xl = 0.360;      % [xk,xl] : [11_ex,12]; to 12 center


x_symmetry = [x0 xa xb xc xd xe xf xg ...
    xh xi xj xk xl];
alpha = zeros(1,length(x_symmetry));
i_symmetry = zeros(length(coils),length(x_symmetry));
%% Calculate Currents at Symmetry points

% Ideally, we want the field gradient in all directions to be constant for
% all positions.  This is not possible given the coil geometry and
% restrictions that the current be positive.

% In every zone of transport (except at the end), we choose so that 
% three coils are active at a time.  This provides a useful constraint,
% because we specify the (B(x),G_z(x), G_x(x)/G_y(x)) at every point.
%
% For boundary matching between each zone, we constrain that one of the
% currents vanishes.  This is convenient, but it requires that we relax
% either G_z(x) or G_x(x)/G_y(x).
%
% The first order calculation will just keep G_z(x) fixed and allow the
% aspect ratio to vary.

G0 = 100;

%%%%%% MOT %%%%%
disp([repmat('%',1,10) ' 0mm ' repmat('%',1,10)]);
kk = 1;
ind = 2;
% Calculate dB_z/dz
dL = 1e-4;
x_this = x_symmetry(kk);
[~,~,Bz1p]=fieldCoil_3D(x_this,0,dL/2,coils(ind).Coil);
[~,~,Bz1n]=fieldCoil_3D(x_this,0,-dL/2,coils(ind).Coil);
Gz = 100*(Bz1p-Bz1n)/dL;
isolve = G0/Gz;
disp([' MOT : ' num2str(round(isolve,1)) ' A']);

alpha(1) = 1; 

i_symmetry(1,1) = isolve;
%%%%%% Transport 3 Coils Sets %%%%%
for kk=2:length(x_symmetry)-1
    x_this = x_symmetry(kk);
    c1 = coils(kk);
    c2 = coils(kk+1);

    disp([repmat('%',1,10) ' ' num2str(x_this*1e3) ' mm ' repmat('%',1,10)]);


    % Calculate B(x)
    [Bx1,By1,Bz1]=fieldCoil_3D(x_this,0,0,c1.Coil);
    [Bx2,By2,Bz2]=fieldCoil_3D(x_this,0,0,c2.Coil);    
    M1 = [Bx1 Bx2;
        By1 By2;
        Bz1 Bz2];
    M1 = M1*1e4; % Convert Tesla to Gauss

    % Calculate dB_x/dx
    dL = 1e-4;
    [Bx1p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,c1.Coil);
    [Bx1n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,c1.Coil);
    Gx1 = 100*(Bx1p-Bx1n)/dL;
    [Bx2p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,c2.Coil);
    [Bx2n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,c2.Coil);
    Gx2 = 100*(Bx2p-Bx2n)/dL;

    % Calculate dB_y/dy
    dL = 1e-4;
    [~,By1p,~]=fieldCoil_3D(x_this,dL/2,0,c1.Coil);
    [~,By1n,~]=fieldCoil_3D(x_this,-dL/2,0,c1.Coil);
    Gy1 = 100*(By1p-By1n)/dL;
    [~,By2p,~]=fieldCoil_3D(x_this,dL/2,0,c2.Coil);
    [~,By2n,~]=fieldCoil_3D(x_this,-dL/2,0,c2.Coil);
    Gy2 = 100*(By2p-By2n)/dL;


    % Calculate dB_z/dz
    dL = 1e-4;
    [~,~,Bz1p]=fieldCoil_3D(x_this,0,dL/2,c1.Coil);
    [~,~,Bz1n]=fieldCoil_3D(x_this,0,-dL/2,c1.Coil);
    Gz1 = 100*(Bz1p-Bz1n)/dL;

    [~,~,Bz2p]=fieldCoil_3D(x_this,0,dL/2,c2.Coil);
    [~,~,Bz2n]=fieldCoil_3D(x_this,0,-dL/2,c2.Coil);
    Gz2 = 100*(Bz2p-Bz2n)/dL;
    M2 = [Gz1 Gz2];

    % Matrix Solve for currents
    b = [0; 0; 0; G0];
    Amat = [M1;M2];

    isolve= mldivide(Amat,b);

    i_symmetry(kk,kk) = isolve(1);
    i_symmetry(kk+1,kk) = isolve(2);


    Gx = Gx1*isolve(1)+Gx2*isolve(2);
    Gy = Gy1*isolve(1)+Gy2*isolve(2);
    Gz = Gz1*isolve(1)+Gz2*isolve(2);

    alpha(kk) = Gy/Gx;
    B = norm(M1*isolve);

    disp([' ' c1.Name ' : ' num2str(round(isolve(1),1)) 'A']);
    disp([' ' c2.Name ' : ' num2str(round(isolve(2),1)) 'A']);
    disp([' (B,Gx,Gy,Gz,Gy/Gx) : ' ...
        '(' num2str(round(B,1)) ' G,' ...
        num2str(round(Gx,1)) ' G/cm,' ...
        num2str(round(Gy,1)) ' G/cm,' ...
        num2str(round(Gz,1)) ',' ...
        num2str(round(Gy/Gx,3)) ')']);
end

%%%%%% Coil 12 %%%%%
x_this = x_symmetry(end);

disp([repmat('%',1,10) ' ' num2str(x_this*1e3) ' mm ' repmat('%',1,10)]);
% Calculate dB_z/dz
dL = 1e-4;
[~,~,Bz1p]=fieldCoil_3D(x_this,0,dL/2,coils(end).Coil);
[~,~,Bz1n]=fieldCoil_3D(x_this,0,-dL/2,coils(end).Coil);
Gz = 100*(Bz1p-Bz1n)/dL;
isolve = G0/Gz;
disp([' Coil 12 : ' num2str(round(isolve,1)) ' A']);
i_symmetry(end,end) = isolve;

alpha(length(coils)) = 1; 

%% Calculate Transport
    figure(20)
    clf

    co=jet(length(coils));


    co = [47 79 79;
        128 0 0;
        0 128 0;
        75, 0, 130;
        255, 140, 0;
255, 255, 0;
0, 255, 0;
0, 255, 255;
0, 0, 255;
255, 0, 255;
238, 232, 170;
100, 149, 237;
255, 105, 180
        ];
co=co/255;
    ps=[];
    legStr={};

%% P,MOT,3
% The coils
ca = coils(1);
cb = coils(2);
cc = coils(3);

% Starting and final aspectr ratio
alpha_i = alpha(1);
alpha_f = alpha(2);

N = 50;
xi = x_symmetry(1);
xf = x_symmetry(2);
X = linspace(xi,xf,N);
A = linspace(alpha_i,alpha_f,N);

i_mat = zeros(3,N);

kk=1;
    for ii = 1:length(X)
        x_this = X(ii);
        alpha_this = A(ii);
        
        % Calculate B(x)
        [By1,Bz1,Bx1]=fieldCoil_3D(0,0,x_this,ca.Coil);
        [Bx2,By2,Bz2]=fieldCoil_3D(x_this,0,0,cb.Coil);   
        [Bx3,By3,Bz3]=fieldCoil_3D(x_this,0,0,cc.Coil);    
     
        M1 = [Bx1 Bx2 Bx3; 
            By1 By2 By3;
            Bz1 Bz2 By3];
        M1 = M1*1e4; % Convert Tesla to Gauss
    
        % Calculate dB_x/dx
        dL = 1e-4;
        [~,~,Bx1p]=fieldCoil_3D(0,0,x_this+dL/2,ca.Coil);
        [~,~,Bx1n]=fieldCoil_3D(0,0,x_this-dL/2,ca.Coil);
        Gx1 = 100*(Bx1p-Bx1n)/dL;
        [Bx2p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,cb.Coil);
        [Bx2n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,cb.Coil);
        Gx2 = 100*(Bx2p-Bx2n)/dL;
        [Bx3p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,cc.Coil);
        [Bx3n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,cc.Coil);
        Gx3 = 100*(Bx3p-Bx3n)/dL;
    
        % Calculate dB_y/dy
        dL = 1e-4;
        [By1p,~,~]=fieldCoil_3D(dL/2,0,x_this,ca.Coil);
        [By1n,~,~]=fieldCoil_3D(-dL/2,0,x_this,ca.Coil);
        Gy1 = 100*(By1p-By1n)/dL;
        [~,By2p,~]=fieldCoil_3D(x_this,dL/2,0,cb.Coil);
        [~,By2n,~]=fieldCoil_3D(x_this,-dL/2,0,cb.Coil);
        Gy2 = 100*(By2p-By2n)/dL;
        [~,By3p,~]=fieldCoil_3D(x_this,dL/2,0,cc.Coil);
        [~,By3n,~]=fieldCoil_3D(x_this,-dL/2,0,cc.Coil);
        Gy3 = 100*(By3p-By3n)/dL;
    
        % Calculate dB_z/dz
        dL = 1e-4;
        [~,Bz1p,~]=fieldCoil_3D(0,dL/2,x_this,ca.Coil);
        [~,Bz1n,~]=fieldCoil_3D(0,-dL/2,x_this,ca.Coil);
        Gz1 = 100*(Bz1p-Bz1n)/dL;
        [~,~,Bz2p]=fieldCoil_3D(x_this,0,dL/2,cb.Coil);
        [~,~,Bz2n]=fieldCoil_3D(x_this,0,-dL/2,cb.Coil);
        Gz2 = 100*(Bz2p-Bz2n)/dL;
        [~,~,Bz3p]=fieldCoil_3D(x_this,0,dL/2,cc.Coil);
        [~,~,Bz3n]=fieldCoil_3D(x_this,0,-dL/2,cc.Coil);
        Gz3 = 100*(Bz3p-Bz3n)/dL; 
    
        M2 = [Gx1 Gx2 Gx3;
            Gy1 Gy2 Gy3;
            Gz1 Gz2 Gz3];
        % M2 = M2*100; % Convert Tesla/m to Gauss/cm
    
    
        %  % Calculate other gradients via Gauss's law
        Gx = -G0/(1+alpha_this);
        Gy = alpha_this*Gx;
         
        b = [0; 0; 0; Gx; Gy; G0];
        Amat = [M1;M2];
        isolve= mldivide(Amat,b);  

        i_mat(1,ii) = isolve(1);
        i_mat(2,ii) = isolve(2);
        i_mat(3,ii) = isolve(3);

        % keyboard
    end



    ps(end+1)=plot(X,i_mat(1,:),'-','color',.7*co(kk,:),'linewidth',2);hold on;
    plot(X,i_mat(2,:),'-','color',.7*co(kk+1,:),'linewidth',2);hold on;
    plot(X,i_mat(3,:),'-','color',.7*co(kk+2,:),'linewidth',2);hold on;
    legStr{end+1}=ca.Name;
    

%%


for kk=2:length(x_symmetry)-2
    % The coils
    ca = coils(kk);
    cb = coils(kk+1);
    cc = coils(kk+2);

        % Starting and final aspectr ratio
    alpha_i = alpha(kk);
    alpha_f = alpha(kk+1);

    N = 50;
    xi = x_symmetry(kk);
    xf = x_symmetry(kk+1);
    X = linspace(xi,xf,N);
    A = linspace(alpha_i,alpha_f,N);

    i_mat = zeros(3,N);

    for ii = 1:length(X)
        x_this = X(ii);
        alpha_this = A(ii);
        
        % Calculate B(x)
        [Bx1,By1,Bz1]=fieldCoil_3D(x_this,0,0,ca.Coil);
        [Bx2,By2,Bz2]=fieldCoil_3D(x_this,0,0,cb.Coil);   
        [Bx3,By3,Bz3]=fieldCoil_3D(x_this,0,0,cc.Coil);    
     
        M1 = [Bx1 Bx2 Bx3; 
            By1 By2 By3;
            Bz1 Bz2 By3];
        M1 = M1*1e4; % Convert Tesla to Gauss
    
        % Calculate dB_x/dx
        dL = 1e-4;
        [Bx1p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,ca.Coil);
        [Bx1n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,ca.Coil);
        Gx1 = 100*(Bx1p-Bx1n)/dL;
        [Bx2p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,cb.Coil);
        [Bx2n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,cb.Coil);
        Gx2 = 100*(Bx2p-Bx2n)/dL;
        [Bx3p,~,~]=fieldCoil_3D(x_this+dL/2,0,0,cc.Coil);
        [Bx3n,~,~]=fieldCoil_3D(x_this-dL/2,0,0,cc.Coil);
        Gx3 = 100*(Bx3p-Bx3n)/dL;
    
        % Calculate dB_y/dy
        dL = 1e-4;
        [~,By1p,~]=fieldCoil_3D(x_this,dL/2,0,ca.Coil);
        [~,By1n,~]=fieldCoil_3D(x_this,-dL/2,0,ca.Coil);
        Gy1 = 100*(By1p-By1n)/dL;
        [~,By2p,~]=fieldCoil_3D(x_this,dL/2,0,cb.Coil);
        [~,By2n,~]=fieldCoil_3D(x_this,-dL/2,0,cb.Coil);
        Gy2 = 100*(By2p-By2n)/dL;
        [~,By3p,~]=fieldCoil_3D(x_this,dL/2,0,cc.Coil);
        [~,By3n,~]=fieldCoil_3D(x_this,-dL/2,0,cc.Coil);
        Gy3 = 100*(By3p-By3n)/dL;
    
        % Calculate dB_z/dz
        dL = 1e-4;
        [~,~,Bz1p]=fieldCoil_3D(x_this,0,dL/2,ca.Coil);
        [~,~,Bz1n]=fieldCoil_3D(x_this,0,-dL/2,ca.Coil);
        Gz1 = 100*(Bz1p-Bz1n)/dL;
        [~,~,Bz2p]=fieldCoil_3D(x_this,0,dL/2,cb.Coil);
        [~,~,Bz2n]=fieldCoil_3D(x_this,0,-dL/2,cb.Coil);
        Gz2 = 100*(Bz2p-Bz2n)/dL;
        [~,~,Bz3p]=fieldCoil_3D(x_this,0,dL/2,cc.Coil);
        [~,~,Bz3n]=fieldCoil_3D(x_this,0,-dL/2,cc.Coil);
        Gz3 = 100*(Bz3p-Bz3n)/dL; 
    
        M2 = [Gx1 Gx2 Gx3;
            Gy1 Gy2 Gy3;
            Gz1 Gz2 Gz3];
        % M2 = M2*100; % Convert Tesla/m to Gauss/cm
    
    
        %  % Calculate other gradients via Gauss's law
        Gx = -G0/(1+alpha_this);
        Gy = alpha_this*Gx;
         
        b = [0; 0; 0; Gx; Gy; G0];
        Amat = [M1;M2];
        isolve= mldivide(Amat,b);  

        i_mat(1,ii) = isolve(1);
        i_mat(2,ii) = isolve(2);
        i_mat(3,ii) = isolve(3);

        % keyboard
    end
% keyboard


    % clf;
        % subplot(1,length(coils),kk);

    ps(end+1)=plot(X,i_mat(1,:),'-','color',.7*co(kk,:),'linewidth',2);hold on;
    plot(X,i_mat(2,:),'-','color',.7*co(kk+1,:),'linewidth',2);hold on;
    plot(X,i_mat(3,:),'-','color',.7*co(kk+2,:),'linewidth',2);hold on;
    legStr{end+1}=ca.Name;
    









end

legend(ps,legStr,'location','eastoutside');