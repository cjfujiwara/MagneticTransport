%% Specify Boundary and Coils
x0 = 0.035;
x1 = 0.068;

alpha0 = 2.5;
alpha1 = 1.707;
G0 = 100;

i1 = 1;
i2 = 2;
i3 = 3;
i4 = 4;
n = 100;

xq = linspace(x0,x1,n);

% Boundary 1
I1 = [I_mat_zone_1(:,end);0];
i1a = I1(1);
i2a = I1(2);
i3a = I1(3);
i4a = I1(4);

% Boundary 2
I3 = I_mat_zone_3(1:4,1);
i1b = I3(1);
i2b = I3(2);
i3b = I3(3);
i4b = I3(4);

%% Construct Set Points and Answers

%Field
B = zeros(1,n);

% Aspect Ratios
A  = linspace(alpha0,alpha1,n);
A = interp1([0.035  0.040 0.064 0.068],[2.5 2.5 1.707 1.707],xq,'pchip');

% Vertical Field Gradient
Gv  = ones(1,n)*G0;

% Horizontal Field Gradient
Gx  = -G0./(1+A);
Gy  = Gx.*A;

% Currents
imat = zeros(4,n);


%% Fields

x_gradient_constraint = Gx;
y_gradient_constraint = Gy;
z_gradient_constraint = Gv;
field_constraint = B;

% Fields
B1 = interp1(X,Bx_all(i1,:),xq,'spline');
B2 = interp1(X,Bx_all(i2,:),xq,'spline');
B3 = interp1(X,Bx_all(i3,:),xq,'spline');
B4 = interp1(X,Bx_all(i4,:),xq,'spline');
field_matrix = [diag(B1(:)) diag(B2(:)) diag(B3(:)) diag(B4(:))];

% Horizontal Gradient
Gx1 = interp1(X,Gx_all(i1,:),xq,'spline');
Gx2 = interp1(X,Gx_all(i2,:),xq,'spline');
Gx3 = interp1(X,Gx_all(i3,:),xq,'spline');
Gx4 = interp1(X,Gx_all(i4,:),xq,'spline');
x_gradient_matrix = [diag(Gx1) diag(Gx2) diag(Gx3) diag(Gx4)];

% Horizontal Gradient
Gy1 = interp1(X,Gy_all(i1,:),xq,'spline');
Gy2 = interp1(X,Gy_all(i2,:),xq,'spline');
Gy3 = interp1(X,Gy_all(i3,:),xq,'spline');
Gy4 = interp1(X,Gy_all(i4,:),xq,'spline');
y_gradient_matrix = [diag(Gy1) diag(Gy2) diag(Gy3) diag(Gy4)];

% Vertical Gradients
Gz1 = interp1(X,Gz_all(i1,:),xq,'spline');
Gz2 = interp1(X,Gz_all(i2,:),xq,'spline');
Gz3 = interp1(X,Gz_all(i3,:),xq,'spline');
Gz4 = interp1(X,Gz_all(i4,:),xq,'spline');
z_gradient_matrix = [diag(Gz1) diag(Gz2) diag(Gz3) diag(Gz4)];


%% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,n*4);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = i1a;               % push
bc_matrix(2,1+n)     = 1; bc_target(2) = i2a;  % MOT
bc_matrix(3,1+2*n)   = 1; bc_target(3) = i3a;               % Coil 3
bc_matrix(4,1+3*n)   = 1; bc_target(4) = i4a;               % Coil 4

bc_matrix(5,n)       = 1; bc_target(5) = i1b;                 % push
bc_matrix(6,2*n)     = 1; bc_target(6) = i2b;                 % MOT
bc_matrix(7,3*n)     = 1; bc_target(7) = i3b;   % Coil 3
bc_matrix(8,4*n)     = 1; bc_target(8) = i4b;   % Coil 4

% % Assemble all constraints
% constraint_matrix = [field_matrix; x_gradient_matrix; y_gradient_matrix; z_gradient_matrix; bc_matrix];
% constraint_vector = [field_constraint(:); x_gradient_constraint(:); y_gradient_constraint(:);z_gradient_constraint(:);bc_target(:)];
% constraint_matrix = sparse(constraint_matrix);

constraint_matrix = [field_matrix; x_gradient_matrix; z_gradient_matrix; bc_matrix];
constraint_vector = [field_constraint(:); x_gradient_constraint(:);z_gradient_constraint(:);bc_target(:)];
constraint_matrix = sparse(constraint_matrix);

%% Initial Guess
i1_guess = linspace(i1a,i1b,n);
i2_guess = linspace(i2a,i2b,n);
i3_guess = linspace(i3a,i3b,n);
i4_guess = linspace(i4a,i4b,n);

init_guess = [i1_guess i2_guess i3_guess i4_guess];

%% Prepare for Optimization

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);
n4 = (1+3*n):(4*n);

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2)+ ...
    sum(diff(curr(n4)).^2);

Imin = 0;
Imax = 150;

lb = Imin*ones(1,numel(init_guess));
% ub = Imax*ones(1,numel(init_guess));
% lb=[];
ub = [];

options.EnableFeasibilityMode = true;
 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,lb,ub,[],options);
 
isolve_1 = x(n1);
isolve_2 = x(n2);
isolve_3 = x(n3);
isolve_4 = x(n4);

% Currents
I_mat_zone_2 = zeros(4,n);

I_mat_zone_2(1,:) = isolve_1;
I_mat_zone_2(2,:) = isolve_2;
I_mat_zone_2(3,:) = isolve_3;
I_mat_zone_2(4,:) = isolve_4;

%% Plot It
hF = figure(12);
clf
hF.Color='w';
hF.Position=[10+350 60 350 850];
subplot(4,1,[1 2 3]);
p1=plot(xq*1e3,I_mat_zone_2(1,:),'color',co(i1,:),'linewidth',2);hold on;
p2=plot(xq*1e3,I_mat_zone_2(2,:),'color',co(i2,:),'linewidth',2);hold on;
p3=plot(xq*1e3,I_mat_zone_2(3,:),'color',co(i3,:),'linewidth',2);hold on;
p4=plot(xq*1e3,I_mat_zone_2(4,:),'color',co(i4,:),'linewidth',2);hold on;

xlabel('position (mm)');
ylabel('current (A)');
set(gca,'fontsize',8,'box','on','linewidth',1);
ylim([0 100]);

xlim([xq(1) xq(end)]*1e3)


legend([p1 p2 p3],legStr([i1 i2 i3]),'location','northwest');

subplot(4,1,4);
c = get(gca,'colororder');
yyaxis left
plot(xq*1e3,Gv,'color',c(1,:),'linewidth',2);
ylim([90 110]);
ylabel('gradient (G/cm)');


yyaxis right
plot(xq*1e3,Gy./Gx,'color',c(2,:),'linewidth',2);hold on;
set(gca,'fontsize',8,'box','on','linewidth',1);
ylabel('aspect ratio');
ylim([.9 3]);
xlim([xq(1) xq(end)]*1e3)
% B0 = isolve_1.*Bx_all(1,:) +  isolve_2.*Bx_all(2,:) + isolve_3.*Bx_all(3,:) + isolve_4.*Bx_all(4,:);
% Gx = isolve_1.*Gx_all(1,:) +  isolve_2.*Gx_all(2,:) + isolve_3.*Gx_all(3,:) + isolve_4.*Gx_all(4,:);
% Gy = isolve_1.*Gy_all(1,:) +  isolve_2.*Gy_all(2,:) + isolve_3.*Gy_all(3,:) + isolve_4.*Gy_all(4,:);
% Gz = isolve_1.*Gz_all(1,:) +  isolve_2.*Gz_all(2,:) + isolve_3.*Gz_all(3,:) + isolve_4.*Gz_all(4,:);
% alpha_me = Gy./Gx;



