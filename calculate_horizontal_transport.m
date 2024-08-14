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

%% Colors

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

%% Construct the coils
G0 = 100;

% Solve vertical transport in the "easy" regime
coils = makeHorizontalCoils;
warning off;

options = optimoptions(@fmincon,'MaxIterations',5e4,'MaxFunctionEvaluations',1e5);

%% Establish Symmetry Points

x0 = 0;0;        % MOT location
xa = 0.0455;0.0455;.044;       % [x0,xa] : [P,MOT,3]; to xa
xb = 0.0680;     % [xa,xb] : [MOT,3,4]; to 3/4 center
xc = 0.1015;     % [xb,xc] : [3,4,5]; to 4/5 center
xd = 0.1310;     % [xc,xd] : [4,5,6]; to 5/6 center
xe = 0.1645;     % [xd,xe] : [5,6,7]; to 6/7 center
xf = 0.1940;     % [xe,xf] : [6,7,8]; to 7/8 center
xg = 0.2275;     % [xf,xg] : [7,8,9]; to 8/9 center
xh = 0.2570;     % [xg,xh] : [8,9,10]; to 9/10 center
xi = 0.2900;0.2959;     % [xh,xi] : [9,10,11]; to 10/11 cetner
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

i_symmetry(2,1) = isolve;
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




    

%% Push Region
% 0 to 68.1 : Push, MOT, Coil 3, Coil 4
N = 50;
X = linspace(x_symmetry(1),x_symmetry(3),N);
dL = 1e-4;

Bx_all = zeros(4,N);
By_all = zeros(4,N);
Bz_all = zeros(4,N);
Gx_all = zeros(4,N);
Gy_all = zeros(4,N);
Gz_all = zeros(4,N);


legStr={};
for kk=1:4
    
    c = coils(kk);
    legStr{kk} = c.Name;

    if kk==1
        [By,Bz,Bx]=fieldCoil_3D(0,0,X,c.Coil);        
        [~,~,Bxp]=fieldCoil_3D(0,0,X+dL/2,c.Coil);
        [~,~,Bxn]=fieldCoil_3D(0,0,X-dL/2,c.Coil);        
        [Byp,~,~]=fieldCoil_3D(dL/2,0,X,c.Coil);
        [Byn,~,~]=fieldCoil_3D(-dL/2,0,X,c.Coil);
        [~,Bzp,~]=fieldCoil_3D(0,dL/2,X,c.Coil);
        [~,Bzn,~]=fieldCoil_3D(0,-dL/2,X,c.Coil);
    else
        [Bx,By,Bz]=fieldCoil_3D(X,0,0,c.Coil);    
        [Bxp,~,~]=fieldCoil_3D(X+dL/2,0,0,c.Coil);
        [Bxn,~,~]=fieldCoil_3D(X-dL/2,0,0,c.Coil);        
        [~,Byp,~]=fieldCoil_3D(X,dL/2,0,c.Coil);
        [~,Byn,~]=fieldCoil_3D(X,-dL/2,0,c.Coil);        
        [~,~,Bzp]=fieldCoil_3D(X,0,dL/2,c.Coil);
        [~,~,Bzn]=fieldCoil_3D(X,0,-dL/2,c.Coil);
    end
    
    % Assign Gradient
    Gx_all(kk,:) = 100*(Bxp-Bxn)/dL;
    Gy_all(kk,:) = 100*(Byp-Byn)/dL;
    Gz_all(kk,:) = 100*(Bzp-Bzn)/dL;
    
    % Assign Field
    Bx_all(kk,:) = 1e4*Bx;
    By_all(kk,:) = 1e4*By;
    Bz_all(kk,:) = 1e4*Bz;    
end

%% Push Region Plot

h1 = figure(10);
clf;
h1.Color='w';

subplot(141);
for kk=1:4
    plot(X*1e3,Bx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('x field (G)');
xlabel('position (mm)');
title('$B_x(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);
legend(legStr,'fontsize',8,'location','best');



subplot(142);
for kk=1:4
    plot(X*1e3,Gx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_x (G/cm)');
xlabel('position (mm)');
title('$\partial_xB_x(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

subplot(143);
for kk=1:4
    plot(X*1e3,Gy_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_y (G/cm)');
xlabel('position (mm)');
title('$\partial_yB_y(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

subplot(144);
for kk=1:4
    plot(X*1e3,Gz_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_z (G/cm)');
xlabel('position (mm)');
title('$\partial_zB_z(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

%% Push Region Optimize

% Target
B0 = 0;    % magnetic field (G)
G0z = G0; % vertical gradient (G/cm)

alpha_i = 1;
alpha_f = 2.5;1.703;
xf = x_symmetry(3);

x2alpha = @(x) (alpha_i-alpha_f)/xf^2*(x-xf).^2+alpha_f;


alpha_vec = x2alpha(X);

alpha_vec = interp1([0 0.02  0.03 .035 .04 .06 xf],[1 2 2.35 2.5 2.5 2 1.703],X,'pchip');

% Field and Gradient Constraints
field_constraint = ones(N,1)*B0;
x_gradient_constraint = -G0z./(1+alpha_vec(:));
y_gradient_constraint = alpha_vec(:).*x_gradient_constraint(:);
z_gradient_constraint = ones(N,1)*G0z;

% field matrix
field_constraint_matrix = [diag(Bx_all(1,:)) diag(Bx_all(2,:)) diag(Bx_all(3,:)) diag(Bx_all(4,:))];
% x-gradient matrix
x_gradient_constraint_matrix = [diag(Gx_all(1,:)) diag(Gx_all(2,:)) diag(Gx_all(3,:)) diag(Gx_all(4,:))];
% y-gradient matrix
y_gradient_constraint_matrix = [diag(Gy_all(1,:)) diag(Gy_all(2,:)) diag(Gy_all(3,:)) diag(Gy_all(4,:))];
% z-gradient matrix
z_gradient_constraint_matrix = [diag(Gz_all(1,:)) diag(Gz_all(2,:)) diag(Gz_all(3,:)) diag(Gz_all(4,:))];


% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,N*3);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = 0;               % push
bc_matrix(2,1+N)     = 1; bc_target(2) = i_symmetry(2,1);  % MOT
bc_matrix(3,1+2*N)   = 1; bc_target(3) = 0;               % Coil 3
bc_matrix(4,1+3*N)   = 1; bc_target(4) = 0;               % Coil 4

bc_matrix(5,N)       = 1; bc_target(5) = 0;                 % push
bc_matrix(6,2*N)     = 1; bc_target(6) = 0;                 % MOT
bc_matrix(7,3*N)     = 1; bc_target(7) = i_symmetry(3,3);   % Coil 3
bc_matrix(8,4*N)     = 1; bc_target(8) = i_symmetry(4,3);   % Coil 4

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; x_gradient_constraint_matrix; y_gradient_constraint_matrix; z_gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; x_gradient_constraint; y_gradient_constraint;z_gradient_constraint;bc_target];
constraint_matrix = sparse(constraint_matrix);



% Initial guess is the simple linear solution

i2a = i_symmetry(2,1);
i3f = i_symmetry(3,3);
i4d = i_symmetry(4,3);

i1_guess = interp1([0 .04 .068],[0 30 0],X,'pchip');
i2_guess = interp1([0 .04 .068],[i2a i2a 0],X,'pchip');
i3_guess = interp1([0 .01 .02 .03 .068],[0 30 60 80 i4d],X,'pchip');
i4_guess = interp1([0 0.02 .04 .068],[0 0 0 i4d],X,'pchip');

% 
% i1_guess = [linspace(0,30,N*.6) linspace(30,0,N*.4)];
% i2_guess = linspace(i_symmetry(2,1),0,N);
% i3_guess = [linspace(0,60,N*.2) linspace(60,i_symmetry(3,3),N*.8)];
% i4_guess = [zeros(1,N*.5) linspace(0,i_symmetry(4,3),N*.5)];
% keyboard

init_guess = [i1_guess i2_guess i3_guess i4_guess];

n1 = (1):(N);
n2 = (1+N):(2*N);
n3 = (1+2*N):(3*N);
n4 = (1+3*N):(4*N);

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2)+ ...
    sum(diff(curr(n4)).^2)+0*sum((curr(n4)<0).*curr(n4).^6);

Imin = 0;
Imax = 150;

lb = Imin*ones(1,numel(init_guess));
% ub = Imax*ones(1,numel(init_guess));
% lb=[];
ub = [];


 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,lb,ub,[],options);
 
isolve_1 = x(n1);
isolve_2 = x(n2);
isolve_3 = x(n3);
isolve_4 = x(n4);

B0 = isolve_1.*Bx_all(1,:) +  isolve_2.*Bx_all(2,:) + isolve_3.*Bx_all(3,:) + isolve_4.*Bx_all(4,:);
Gx = isolve_1.*Gx_all(1,:) +  isolve_2.*Gx_all(2,:) + isolve_3.*Gx_all(3,:) + isolve_4.*Gx_all(4,:);
Gy = isolve_1.*Gy_all(1,:) +  isolve_2.*Gy_all(2,:) + isolve_3.*Gy_all(3,:) + isolve_4.*Gy_all(4,:);
Gz = isolve_1.*Gz_all(1,:) +  isolve_2.*Gz_all(2,:) + isolve_3.*Gz_all(3,:) + isolve_4.*Gz_all(4,:);
alpha_me = Gy./Gx;


% G_0a = i1_0a.*g1 + i2_0a.*g2 + i3_0a.*g3;
% B_0a = i1_0a.*b1 + i2_0a.*b2 + i3_0a.*g3;

%% Plot Solved Values

h2 = figure(11);
clf;
h2.Color='w';

subplot(3,1,[1 2]);
plot(X*1e3,isolve_1,'color',co(1,:),'linewidth',2);hold on
plot(X*1e3,isolve_2,'color',co(2,:),'linewidth',2);
plot(X*1e3,isolve_3,'color',co(3,:),'linewidth',2);
plot(X*1e3,isolve_4,'color',co(4,:),'linewidth',2);
xlabel('position (mm)');
ylabel('current (A)');
xlim([x_symmetry(1) x_symmetry(3)]*1e3);
legend(legStr,'fontsize',8,'location','best');
set(gca,'fontsize',10,'box','on','linewidth',1,'xaxislocation','top');
ylim([0 100]);

yyaxis right
plot(X*1e3,B0);
ylim([-1 1]);
yyaxis left


subplot(3,1,3);
c2=get(gca,'colororder');
yyaxis left
plot(X*1e3,Gz,'color',c2(1,:),'linewidth',2);hold on
ylabel('z gradient (G/cm)');
xlim([x_symmetry(1) x_symmetry(3)]*1e3);
ylim([99 101]);

yyaxis right
plot(X*1e3,alpha_me,'color',c2(2,:),'linewidth',2);hold on
ylabel('horizontal aspect ratio');
ylim([.9 3]);
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('position (mm)');


% plot(X,isolve_2,'color',co(2,:),'linewidth',2);
% plot(X,isolve_3,'color',co(3,:),'linewidth',2);
% plot(X,isolve_4,'color',co(4,:),'linewidth',2);
%% "Normal" Zones
i1 = 2;
i2 = 9;
legStr = {};

    hF_normal = figure(20);
    hF_normal.Color='w';
    clf
ps=[];
for kk=i1:i2
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
    end

    ps(end+1)=plot(1e3*X,i_mat(1,:),'-','color',.7*co(kk,:),'linewidth',2);hold on;
    ps(end+1)=plot(1e3*X,i_mat(2,:),'-','color',.7*co(kk+1,:),'linewidth',2);hold on;
    ps(end+1)=plot(1e3*X,i_mat(3,:),'-','color',.7*co(kk+2,:),'linewidth',2);hold on;
    legStr{end+1}=ca.Name;
    legStr{end+1}=cb.Name;
    legStr{end+1}=cc.Name;        
end

[~,sinds]=unique(legStr,'stable');
legend(ps(sinds),legStr(sinds),'location','eastoutside');
xlim([x_symmetry(i1) x_symmetry(i2+1)]*1e3);
xlabel('position (mm)');
ylabel('current (A)');
ylim([0 70]);
title('Normal Transport');


    

%% HV Transfer Region
% 290 to 360 : Coil 10, Coil 11, Coil 11 Extra, Coil 12
N = 26;
X = linspace(x_symmetry(10),x_symmetry(end),N);
dL = 1e-4;

Bx_all = zeros(4,N);
By_all = zeros(4,N);
Bz_all = zeros(4,N);
Gx_all = zeros(4,N);
Gy_all = zeros(4,N);
Gz_all = zeros(4,N);

inds = [10 11 12 13];
legStr={};
for kk=1:4
    
    c = coils(inds(kk));
    legStr{kk} = c.Name;

    [Bx,By,Bz]=fieldCoil_3D(X,0,0,c.Coil);    
    [Bxp,~,~]=fieldCoil_3D(X+dL/2,0,0,c.Coil);
    [Bxn,~,~]=fieldCoil_3D(X-dL/2,0,0,c.Coil);        
    [~,Byp,~]=fieldCoil_3D(X,dL/2,0,c.Coil);
    [~,Byn,~]=fieldCoil_3D(X,-dL/2,0,c.Coil);        
    [~,~,Bzp]=fieldCoil_3D(X,0,dL/2,c.Coil);
    [~,~,Bzn]=fieldCoil_3D(X,0,-dL/2,c.Coil);
    
    % Assign Gradient
    Gx_all(kk,:) = 100*(Bxp-Bxn)/dL;
    Gy_all(kk,:) = 100*(Byp-Byn)/dL;
    Gz_all(kk,:) = 100*(Bzp-Bzn)/dL;
    
    % Assign Field
    Bx_all(kk,:) = 1e4*Bx;
    By_all(kk,:) = 1e4*By;
    Bz_all(kk,:) = 1e4*Bz;    
end

%% HV Region Plot

h1 = figure(30);
clf;
h1.Color='w';

subplot(141);
for kk=1:4
    plot(X*1e3,Bx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('x field (G)');
xlabel('position (mm)');
title('$B_x(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);
legend(legStr,'fontsize',8,'location','best');



subplot(142);
for kk=1:4
    plot(X*1e3,Gx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_x (G/cm)');
xlabel('position (mm)');
title('$\partial_xB_x(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

subplot(143);
for kk=1:4
    plot(X*1e3,Gy_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_y (G/cm)');
xlabel('position (mm)');
title('$\partial_yB_y(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

subplot(144);
for kk=1:4
    plot(X*1e3,Gz_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_z (G/cm)');
xlabel('position (mm)');
title('$\partial_zB_z(x)$','interpreter','latex');
set(gca,'FontSize',12,'box','on','linewidth',1);

%% HV Region Optimize

% Target
B0 = 0;    % magnetic field (G)
G0z = G0; % vertical gradient (G/cm)

alpha_i = 1.703;
alpha_f = 1;


alpha_vec = interp1([.290 .340 .360],[1.703 1.703 1],X,'pchip');

% Field and Gradient Constraints
field_constraint = ones(N,1)*B0;
x_gradient_constraint = -G0z./(1+alpha_vec(:));
y_gradient_constraint = alpha_vec(:).*x_gradient_constraint(:);
z_gradient_constraint = ones(N,1)*G0z;

% field matrix
field_constraint_matrix = [diag(Bx_all(1,:)) diag(Bx_all(2,:)) diag(Bx_all(3,:)) diag(Bx_all(4,:))];
% x-gradient matrix
x_gradient_constraint_matrix = [diag(Gx_all(1,:)) diag(Gx_all(2,:)) diag(Gx_all(3,:)) diag(Gx_all(4,:))];
% y-gradient matrix
y_gradient_constraint_matrix = [diag(Gy_all(1,:)) diag(Gy_all(2,:)) diag(Gy_all(3,:)) diag(Gy_all(4,:))];
% z-gradient matrix
z_gradient_constraint_matrix = [diag(Gz_all(1,:)) diag(Gz_all(2,:)) diag(Gz_all(3,:)) diag(Gz_all(4,:))];


% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,N*3);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = i_symmetry(10,10);               % coil 10
bc_matrix(2,1+N)     = 1; bc_target(2) = i_symmetry(11,10);  % coil 11
bc_matrix(3,1+2*N)   = 1; bc_target(3) = 0;               % Coil 11 extra
bc_matrix(4,1+3*N)   = 1; bc_target(4) = 0;               % Coil 12

bc_matrix(5,N)       = 1; bc_target(5) = 0;                  % coil 10
bc_matrix(6,2*N)     = 1; bc_target(6) = 0;                  % coil 11
bc_matrix(7,3*N)     = 1; bc_target(7) = 0;                  % coil 11 extra
bc_matrix(8,4*N)     = 1; bc_target(8) = i_symmetry(end,end);% coil 12

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; x_gradient_constraint_matrix; y_gradient_constraint_matrix; z_gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; x_gradient_constraint; y_gradient_constraint;z_gradient_constraint;bc_target];
constraint_matrix = sparse(constraint_matrix);



% Initial guess is the simple linear solution
i1a = i_symmetry(10,10);
i2a = i_symmetry(11,10);

i2b = i_symmetry(11,11);
i3b = i_symmetry(12,11);

i3c = i_symmetry(12,12);
i4c = i_symmetry(13,12);


i4d = i_symmetry(end,end);

i1_guess = interp1(x_symmetry(10:13),[i1a 0 0 0],X,'linear');
i2_guess = interp1(x_symmetry(10:13),[i2a i2b 0 0],X,'linear');
i3_guess = interp1(x_symmetry(10:13),[0  i3b i3c 0],X,'linear');
i4_guess = interp1(x_symmetry(10:13),[0 0 i4c   i4d],X,'linear');

% 
% i1_guess = [linspace(0,30,N*.6) linspace(30,0,N*.4)];
% i2_guess = linspace(i_symmetry(2,1),0,N);
% i3_guess = [linspace(0,60,N*.2) linspace(60,i_symmetry(3,3),N*.8)];
% i4_guess = [zeros(1,N*.5) linspace(0,i_symmetry(4,3),N*.5)];
% keyboard

init_guess = [i1_guess i2_guess i3_guess i4_guess];

n1 = (1):(N);
n2 = (1+N):(2*N);
n3 = (1+2*N):(3*N);
n4 = (1+3*N):(4*N);

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2)+ ...
    sum(diff(curr(n4)).^2)+0*sum((curr(n4)<0).*curr(n4).^6);

Imin = 0;
Imax = 150;

lb = Imin*ones(1,numel(init_guess));
% ub = Imax*ones(1,numel(init_guess));
% lb=[];
ub = [];
options = optimoptions(@fmincon,'MaxIterations',5e4,'MaxFunctionEvaluations',1e6);

disp('optimizing current traces ....');
 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,lb,ub,[],options);
 disp('done');
 
isolve_1 = x(n1);
isolve_2 = x(n2);
isolve_3 = x(n3);
isolve_4 = x(n4);

B0 = isolve_1.*Bx_all(1,:) +  isolve_2.*Bx_all(2,:) + isolve_3.*Bx_all(3,:) + isolve_4.*Bx_all(4,:);
Gx = isolve_1.*Gx_all(1,:) +  isolve_2.*Gx_all(2,:) + isolve_3.*Gx_all(3,:) + isolve_4.*Gx_all(4,:);
Gy = isolve_1.*Gy_all(1,:) +  isolve_2.*Gy_all(2,:) + isolve_3.*Gy_all(3,:) + isolve_4.*Gy_all(4,:);
Gz = isolve_1.*Gz_all(1,:) +  isolve_2.*Gz_all(2,:) + isolve_3.*Gz_all(3,:) + isolve_4.*Gz_all(4,:);
alpha_me = Gy./Gx;


% G_0a = i1_0a.*g1 + i2_0a.*g2 + i3_0a.*g3;
% B_0a = i1_0a.*b1 + i2_0a.*b2 + i3_0a.*g3;

%% Plot Solved Values

h2 = figure(31);
clf;
h2.Color='w';

subplot(3,1,[1 2]);
plot(X*1e3,isolve_1,'color',co(1,:),'linewidth',2);hold on
plot(X*1e3,isolve_2,'color',co(2,:),'linewidth',2);
plot(X*1e3,isolve_3,'color',co(3,:),'linewidth',2);
plot(X*1e3,isolve_4,'color',co(4,:),'linewidth',2);
xlabel('position (mm)');
ylabel('current (A)');
xlim([x_symmetry(10) x_symmetry(13)]*1e3);
legend(legStr,'fontsize',8,'location','best');
set(gca,'fontsize',10,'box','on','linewidth',1,'xaxislocation','top');
ylim([0 100]);

yyaxis right
plot(X*1e3,B0);
ylim([-1 1]);
yyaxis left


subplot(3,1,3);
c2=get(gca,'colororder');
yyaxis left
plot(X*1e3,Gz,'color',c2(1,:),'linewidth',2);hold on
ylabel('z gradient (G/cm)');
xlim([x_symmetry(10) x_symmetry(13)]*1e3);
ylim([99 101]);

yyaxis right
plot(X*1e3,alpha_me,'color',c2(2,:),'linewidth',2);hold on
ylabel('horizontal aspect ratio');
ylim([.9 3]);
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('position (mm)');


% plot(X,isolve_2,'color',co(2,:),'linewidth',2);
% plot(X,isolve_3,'color',co(3,:),'linewidth',2);
% plot(X,isolve_4,'color',co(4,:),'linewidth',2);