%% Specify Boundary and Coils
x0 = 0.330;
x1 = 0.360;

alpha0 = 1.707;
alpha1 = 1;
G0 = 100;

i1 = 11;
i2 = 12;
i3 = 13;
n = 100;

xq = linspace(x0,x1,n);

% Boundary 1
I1 = I_mat_zone_4(:,end);
i1a = I1(2);
i2a = I1(3);
i3a = I1(4);


% Boundary 2
i1b = 0;
i2b = 0;
i3b =  18.8887;

%% Construct Set Points and Answers

%Field
B = zeros(1,n);

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
field_matrix = [diag(B1(:)) diag(B2(:)) diag(B3(:))];

% Vertical Gradients
Gz1 = interp1(X,Gz_all(i1,:),xq,'spline');
Gz2 = interp1(X,Gz_all(i2,:),xq,'spline');
Gz3 = interp1(X,Gz_all(i3,:),xq,'spline');
z_gradient_matrix = [diag(Gz1) diag(Gz2) diag(Gz3)];


%% Construct the boundary condition constraint matrix
bc_matrix = zeros(6,n*3);
bc_target = zeros(6,1);
bc_matrix(1,1)       = 1; bc_target(1) = i1a;               % push
bc_matrix(2,1+n)     = 1; bc_target(2) = i2a;  % MOT
bc_matrix(3,1+2*n)   = 1; bc_target(3) = i3a;               % Coil 3

bc_matrix(4,n)       = 1; bc_target(4) = i1b;                 % push
bc_matrix(5,2*n)     = 1; bc_target(5) = i2b;                 % MOT
bc_matrix(6,3*n)     = 1; bc_target(6) = i3b;   % Coil 3

% % Assemble all constraints
% constraint_matrix = [field_matrix; x_gradient_matrix; y_gradient_matrix; z_gradient_matrix; bc_matrix];
% constraint_vector = [field_constraint(:); x_gradient_constraint(:); y_gradient_constraint(:);z_gradient_constraint(:);bc_target(:)];
% constraint_matrix = sparse(constraint_matrix);

constraint_matrix = [field_matrix; z_gradient_matrix; bc_matrix];
constraint_vector = [field_constraint(:);z_gradient_constraint(:);bc_target(:)];
constraint_matrix = sparse(constraint_matrix);

%% Initial Guess
i1_guess = linspace(i1a,i1b,n);
i2_guess = linspace(i2a,i2b,n);
i3_guess = linspace(i3a,i3b,n);

init_guess = [i1_guess i2_guess i3_guess];

%% Prepare for Optimization

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2);

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

B = isolve_1.*B1+isolve_2.*B2+isolve_3.*B3;
Gx = isolve_1.*Gx1+isolve_2.*Gx2+isolve_3.*Gx3;
Gy = isolve_1.*Gy1+isolve_2.*Gy2+isolve_3.*Gy3;
Gz = isolve_1.*Gz1+isolve_2.*Gz2+isolve_3.*Gz3;





% Currents
I_mat_zone_5b = zeros(4,n);

I_mat_zone_5b(1,:) = isolve_1;
I_mat_zone_5b(2,:) = isolve_2;
I_mat_zone_5b(3,:) = isolve_3;

%% Plot It
hF = figure(16);
clf
hF.Color='w';
hF.Position=[10+350 60 350 850];
subplot(4,1,[1 2 3]);
p1=plot(xq*1e3,I_mat_zone_5b(1,:),'color',co(i1,:),'linewidth',2);hold on;
p2=plot(xq*1e3,I_mat_zone_5b(2,:),'color',co(i2,:),'linewidth',2);hold on;
p3=plot(xq*1e3,I_mat_zone_5b(3,:),'color',co(i3,:),'linewidth',2);hold on;

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



