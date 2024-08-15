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

%% Calculate Fields and Gradients
fprintf('Calculating magnetic field profile for all coils ...')
N = 10000;
X = linspace(0,0.360,N);
dL = 1e-4;

Bx_all = zeros(length(coils),N);
By_all = zeros(length(coils),N);
Bz_all = zeros(length(coils),N);
Gx_all = zeros(length(coils),N);
Gy_all = zeros(length(coils),N);
Gz_all = zeros(length(coils),N);
legStr={};
for kk=1:length(coils)   
    c = coils(kk);
    legStr{kk}=c.Name;

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

disp(' done!');

%% Plot all Coils

h1 = figure(10);
clf;
h1.Color='w';
h1.Position = [10 60 850 850];

ax1=subplot(411);
for kk=1:length(coils)
    plot(X*1e3,Bx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('x field (G)');
xlabel('position (mm)');
title('$B_x(x)$','interpreter','latex');
set(gca,'FontSize',10,'box','on','linewidth',1);
legend(legStr,'fontsize',6,'location','eastoutside');
xlim([X(1) X(end)]*1e3);


ax2=subplot(412);
for kk=1:length(coils)
    plot(X*1e3,Gx_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_x (G/cm)');
xlabel('position (mm)');
title('$\partial_xB_x(x)$','interpreter','latex');
set(gca,'FontSize',10,'box','on','linewidth',1);
legend(legStr,'fontsize',6,'location','eastoutside');
xlim([X(1) X(end)]*1e3);

ax3=subplot(413);
for kk=1:length(coils)
    plot(X*1e3,Gy_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_y (G/cm)');
xlabel('position (mm)');
title('$\partial_yB_y(x)$','interpreter','latex');
set(gca,'FontSize',10,'box','on','linewidth',1);
legend(legStr,'fontsize',6,'location','eastoutside');
xlim([X(1) X(end)]*1e3);

ax4=subplot(414);
for kk=1:length(coils)
    plot(X*1e3,Gz_all(kk,:),'-','color',co(kk,:),'linewidth',2); 
    hold on;
end
ylabel('G_z (G/cm)');
xlabel('position (mm)');
title('$\partial_zB_z(x)$','interpreter','latex');
set(gca,'FontSize',10,'box','on','linewidth',1);
legend(legStr,'fontsize',6,'location','eastoutside');
xlim([X(1) X(end)]*1e3);

linkaxes([ax1 ax2 ax3 ax4],'x');

%% Zone 1 Calculation

calculate_horizontal_transport_zone_1
calculate_horizontal_transport_zone_3
calculate_horizontal_transport_zone_5

calculate_horizontal_transport_zone_2
calculate_horizontal_transport_zone_4
 
% calculate_horizontal_transport_zone_5b
