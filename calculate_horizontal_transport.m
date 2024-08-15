%% Introduction
% calculate_horizontal_transport.m
%
% Author : CJ Fujiwara
%
% This code caluclates the current traces for horizontal transport which 
% smoothly transforms a magnetic trap formed at the center fo the MOT
% chamber to one that is formed at the center of Coil 12.
%
% The calculation is divided into four zones each of which are optimized
% differently.
%
% zone1 : MOT to  MOT edge; uses [Push,MOT,Coil 3]
% zone2 : MOT edge to Coil 3 edge; uses [Push,MOT,Coil 3,Coil 4]
% zone3 : Coil 3 edge to Coil 11 edge; uses Coils 3 to 11
% zone4 : Coil 11 edge to Coil 12 edge; [Coil 10, Coil 11, Coil 11 extra, Coil 12]
% zone5 : Coil 12 edge to Coil 12; [Coil 11, Coil 11 extra, Coil 12]
%
% In Zone 1,3, and 5, only three coils are active at a time, in which an
% exact solution is specified by three constraints (field zero, vertical
% gradient, horizontal aspect ratio).
%
% In Zone 2 and 4, four coils are active.  These regions uses four coils
% because the fields produced are relatively weak because of vacuum chamber
% constraints.  The problem is overconstrained and so a optimizer is
% utilized which minimizes the current curvature integral( (dI/dx)^2*dx) 
% which reduces the 'kinkiness' in the current traces.
%
%
% makeHorizontalCoils - creates a structure array which contains
% descriptions for every single coil.
%
% Notes :
% - negative currents are disallowed due to restrictions on the
% electronics. If negative currents are possible, the currents traces will
% get much better and have a better aspect ratio.
%
% References :
% Michael Yee Thesis
% Michael Yee Data folder
% Dave McKay Thesis
% Dave McKay Data folder

%% Settings
% Target vertical field gradient [G/cm]
G0 = 100; 

doSave = 1;

% Construct the coil object
coils = makeHorizontalCoils;
warning off;

% Options for current trace optimizer
options = optimoptions(@fmincon,'MaxIterations',5e4,...
    'MaxFunctionEvaluations',1e5);
%% Add to Path
% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

chdir(curpath);

a = fileparts(curpath);
addpath(a);addpath(genpath(a));
%% Colors
% Unique color vector for 13 coils
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
    255, 105, 180]/255;

%% Calculate Fields and Gradients
% Calculate all magnetic fields
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
% Plot field profile for every single coil
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

%% Calculate Transport in Each Zone
% This is main part of the code which is divided in five parts. Caluclate
% the current traces for each zone.

calculate_horizontal_transport_zone_1;
calculate_horizontal_transport_zone_3;
calculate_horizontal_transport_zone_5;

% Require optimization
calculate_horizontal_transport_zone_2;
calculate_horizontal_transport_zone_4;
 
% calculate_horizontal_transport_zone_5b

%% Assemble Calculations and Resample
% Assemble the results of each calculation and convert into a single
% object

xq_p = [xq1 xq2 xq3 xq4 xq5];
I_mat_zone_1_p = [I_mat_zone_1; zeros(10,size(I_mat_zone_1,2))];
I_mat_zone_2_p = [I_mat_zone_2; zeros(9,size(I_mat_zone_2,2))];

I_mat_zone_4_p = [zeros(9,size(I_mat_zone_4,2)); I_mat_zone_4];
I_mat_zone_5_p = [zeros(10,size(I_mat_zone_5,2));I_mat_zone_5];

I_mat_all = [I_mat_zone_1_p I_mat_zone_2_p I_mat_zone_3 I_mat_zone_4_p I_mat_zone_5_p ];

X_out = 0:1e-4:0.360; 
I_out = zeros(13,length(X_out));

for rr=1:13
    x = xq_p;
    y = I_mat_all(rr,:);
    [xu,inds] = unique(x);
    yu = y(inds);
    I_out(rr,:)=interp1(xu,yu,X_out,'linear');
end

% Recalculate the fields and gradients
legStr={};dL = 1e-4;

Bx0 = zeros(length(coils),numel(X_out));
By0 = zeros(length(coils),numel(X_out));
Bz0 = zeros(length(coils),numel(X_out));
Gx0 = zeros(length(coils),numel(X_out));
Gy0 = zeros(length(coils),numel(X_out));
Gz0 = zeros(length(coils),numel(X_out));

for kk=1:length(coils)   
    c = coils(kk);
    legStr{kk}=c.Name;

    if kk==1
        [By,Bz,Bx]=fieldCoil_3D(0,0,X_out,c.Coil);        
        [~,~,Bxp]=fieldCoil_3D(0,0,X_out+dL/2,c.Coil);
        [~,~,Bxn]=fieldCoil_3D(0,0,X_out-dL/2,c.Coil);        
        [Byp,~,~]=fieldCoil_3D(dL/2,0,X_out,c.Coil);
        [Byn,~,~]=fieldCoil_3D(-dL/2,0,X_out,c.Coil);
        [~,Bzp,~]=fieldCoil_3D(0,dL/2,X_out,c.Coil);
        [~,Bzn,~]=fieldCoil_3D(0,-dL/2,X_out,c.Coil);
    else
        [Bx,By,Bz]=fieldCoil_3D(X_out,0,0,c.Coil);    
        [Bxp,~,~]=fieldCoil_3D(X_out+dL/2,0,0,c.Coil);
        [Bxn,~,~]=fieldCoil_3D(X_out-dL/2,0,0,c.Coil);        
        [~,Byp,~]=fieldCoil_3D(X_out,dL/2,0,c.Coil);
        [~,Byn,~]=fieldCoil_3D(X_out,-dL/2,0,c.Coil);        
        [~,~,Bzp]=fieldCoil_3D(X_out,0,dL/2,c.Coil);
        [~,~,Bzn]=fieldCoil_3D(X_out,0,-dL/2,c.Coil);
    end    
    % Assign Gradient
    Gx0(kk,:) = 100*(Bxp-Bxn)/dL;
    Gy0(kk,:) = 100*(Byp-Byn)/dL;
    Gz0(kk,:) = 100*(Bzp-Bzn)/dL;    
    % Assign Field
    Bx0(kk,:) = 1e4*Bx;
    By0(kk,:) = 1e4*By;
    Bz0(kk,:) = 1e4*Bz;    
end

Gx_out = sum(Gx0.*I_out,1);
Gy_out = sum(Gy0.*I_out,1);
Gz_out = sum(Gz0.*I_out,1);
Bx_out = sum(Bx0.*I_out,1);
By_out = sum(By0.*I_out,1);
Bz_out = sum(Bz0.*I_out,1);

A_out = Gy_out./Gx_out;


%% Plot it

ff=figure(100);
ff.Color='w';
ff.Position=[50 50 900 600];
clf
c= get(gca,'colororder');
subplot(4,1,[1 2 3]);
legStr={};
for kk=1:13
    legStr{kk}=coils(kk).Name;

plot(X_out*1e3,I_out(kk,:),'color',co(kk,:),'linewidth',2); hold on
end
xlim([X_out(1) X_out(end)]*1e3)
ylabel('current (A)');
legend(legStr,'fontsize',6,'location','north','numcolumns',5);
xlabel('position (mm)');
set(gca,'box','on','linewidth',1,'fontsize',10,'xaxislocation','top');
ylim([-1 100]);

subplot(4,1,[4]);
yyaxis left
plot(X_out*1e3,Gz_out,'color',c(1,:),'linewidth',2);hold on
ylim([90 110]);
ylabel('gradient (G/cm)');
yyaxis right
plot(X_out*1e3,A_out,'color',c(2,:),'linewidth',2)
ylim([.9 3]);
ylabel('aspect ratio');
xlim([X_out(1) X_out(end)]*1e3)
xlabel('position (mm)');
set(gca,'box','on','linewidth',1,'fontsize',10);

%% Export the Output

if doSave
    % Add all subdirectories for this m file
    curpath = fileparts(mfilename('fullpath'));
    outdir = fullfile(curpath,'horizontal_output');       
    saveas(h1,fullfile(outdir,'horizontal_field_profile.png'));
    saveas(h1,fullfile(outdir,'horizontal_field_profile.fig'));
    saveas(ff,fullfile(outdir,'horizontal_current.png'));
    saveas(ff,fullfile(outdir,'horizontal_current.fig'));

    save(fullfile(outdir,'horizontal_current'),'X_out','I_out','Gx_out','Gy_out','Gz_out','Bx_out','By_out','Bz_out','A_out');
end
