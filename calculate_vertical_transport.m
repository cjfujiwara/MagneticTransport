%% Add to Path

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
coils = makeVerticalCoils;
warning off;

options = optimoptions(@fmincon,'MaxIterations',5e4,'MaxFunctionEvaluations',1e5);

%% Establish Symmetry Points
% During vertical transport there are 6 points of high symmetry. 

% The first class of high symmetry points are at the beginning and end 
% of transport where the trap is formed at the center of top two (15 + 16) 
% or bottom two coils (12a + 12b). This is important as these are starting
% stage from horizontal transport and the magnetic evaporation. Although
% the magnetic evaporation can occur at an offset position to avoid hitting
% the sapphire window
%
% The second class of high points are at the centers of each of the
% interior coils (12b, 13, 14, 15). These points are of particular interest
% because at those heights, that particular coils adds no gradient.  So it
% is simple to have its neighboring coils provide the requisite gradient.
% For example, at the height of Coil 13, only Coil 12b and Coil 14 should
% be providing a magnetic gradient.  These high symmetry points create the
% natural "cell" for calculating the transport.
%
% Because the problem is over contstrained, we shall constrain the problem
% by minizing the current curvature which makes the current regulation the
% smoothest.
%
% To minimize the curvature we desire to minimze the function
%
% func({I_i}(z)) = sum_i ( integral_za_to_zb ( dI_i/dz)^2 dz)
% Where I_i(z) is the current for a given coil and za and zb are the edges
% of a unit cell.  This is a very simple minimization problem whose
% solution is just a line (if all coils and intervals are symmetric) 
% 
% In the presnece of realistic assymetries, the boundary conditions modify
% the solution to this minimization problem
%
% The constraint is a simple linear constraint on all I_i(z) such that the
% field gradient is some target value G0 and the B field is some target
% value B0 = 0 Gauss.
%
% The final constraint are the boundary conditions which we set by the high
% symmetry points which are physically motivated.

% Find center position of each coil
z_centers=zeros(6,1);
for ii=1:length(coils)
    z_centers(ii) = coils(ii).zbot+coils(ii).Height/2;
end

% Idenfity High Symmetry Points : calculate
z_init = 0;
za = z_centers(2);
zb = z_centers(3);
zc = z_centers(4);
zd = z_centers(5);
z_final = (z_centers(5)+z_centers(6))*.5;

% Idenfity High Symmetry Points : auto specify
z_init = 0; 
za = 0.023; 
zb = 0.068;
zc = 0.102;
zd = 0.148;
z_final = (z_centers(5)+z_centers(6))*.5;

z_symmetry = [z_init za zb zc zd z_final];
%% Initialize data vectors
N = 1e4; % Number of positions to evalulate
dL  = 1e-4; % Distance separation for calculating the gradient
Z = linspace(z_init,z_final,N);

B0 = 0;     % Target field in Gauss
G0 = 80;     % Target gradient in Gauss/cm

% Vectors for magnetic fields and gradients
B = zeros(N,6);
G = zeros(N,6);
%% Calculate the field and gradient at all points in space for all coils

for kk=1:length(coils)
    C = coils(kk).Coil;    
    [~,~,Bz]=fieldCoil_3D(0,0,Z,C);
    [~,~,Bzp]=fieldCoil_3D(0,0,Z+dL,C);
    [~,~,Bzn]=fieldCoil_3D(0,0,Z-dL,C);
    Gz = (Bzp-Bzn)/(2*dL);
    B(:,kk) = Bz;
    G(:,kk) = Gz;
end

% Convert to Gauss and Gauss/cm
B = B*1e4;
G = G*1e2;

%% Calculate Field and Gradient at symmetry points
B_bc = zeros(6,6);
G_bc = zeros(6,6);

I_bc = zeros(6,6);

for kk=1:length(coils)
    B_bc(:,kk) = interp1(Z,B(:,kk)',z_symmetry,'spline');
    G_bc(:,kk) = interp1(Z,G(:,kk)',z_symmetry,'spline');
end

%% Calculate the set currents at the high symmetry points

% at z_init only want Coil 1 and Coil 2
A = [B_bc(1,1) B_bc(1,2);
    G_bc(1,1) G_bc(1,2)];
i0=A\[B0; G0];

% at za only want Coil 1 and Coil 3
A = [B_bc(2,1) B_bc(2,3);
    G_bc(2,1) G_bc(2,3)];
ia=A\[B0; G0];

% at zb only want Coil 2 and Coil 4
A = [B_bc(3,2) B_bc(3,4);
    G_bc(3,2) G_bc(3,4)];
ib=A\[B0; G0];

% at zc only want Coil 3 and Coil 5
A = [B_bc(4,3) B_bc(4,5);
    G_bc(4,3) G_bc(4,5)];
ic=A\[B0; G0];

% at zd only want Coil 4 and Coil 6
A = [B_bc(5,4) B_bc(5,6);
    G_bc(5,4) G_bc(5,6)];
id=A\[B0; G0];

% at zfinal only want Coil 5 and Coil 6
A = [B_bc(6,5) B_bc(6,6);
    G_bc(6,5) G_bc(6,6)];
ifinal=A\[B0; G0];


%% Plot Field And Gradient
figure(1);
clf
co=get(gca,'colororder');

subplot(211)
for kk=1:length(coils)
    plot(Z*1e3,B(:,kk),'-','linewidth',1,'color',co(kk,:));
    hold on
end
xlabel('position (mm)')
ylabel('B(z) @ 1A (G)')
legend({'12a','12b','13','14','15','16'},'orientation','horizontal','location','northwest')
xlim([min(Z) max(Z)]*1e3);

subplot(212)
for kk=1:length(coils)
    plot(Z*1e3,G(:,kk),'-','linewidth',1,'color',co(kk,:));
    hold on
end
xlabel('position (mm)')
ylabel('dB/dz @ 1A (G/cm)')
legend({'12a','12b','13','14','15','16'},'orientation','horizontal','location','southwest')
xlim([min(Z) max(Z)]*1e3);

%% Zone init to a calculation
n = 100;                % Points in this zone to evaluate
z_0a = linspace(z_init,za,n);  % Vector of points to calculate

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);

% Recalculate the field at this specific mesh because it is numerically
% easier to solve the problem on a mesh that is equally spaced between the
% high symmetry points.

% Find the Bfield at these points
b1 = interp1(Z,B(:,1),z_0a,'spline');
b2 = interp1(Z,B(:,2),z_0a,'spline');
b3 = interp1(Z,B(:,3),z_0a,'spline');

% Find the gradient at these points
g1 = interp1(Z,G(:,1),z_0a,'spline');
g2 = interp1(Z,G(:,2),z_0a,'spline');
g3 = interp1(Z,G(:,3),z_0a,'spline');

% Construct the field constraint matrix
field_constraint_matrix = [diag(b1) diag(b2) diag(b3)];
field_constraint = ones(n,1)*B0;

% Construct the gradient constraint matrix
gradient_constraint_matrix = [diag(g1) diag(g2) diag(g3)];
gradient_constraint = ones(n,1)*G0;

% Construct the boundary condition constraint matrix
bc_matrix = zeros(6,n*3);
bc_target = zeros(6,1);
bc_matrix(1,1)       = 1; bc_target(1) = i0(1);     % I1(za) = calculated
bc_matrix(2,1+n)     = 1; bc_target(2) = i0(2);     % I2(za) = calculated
bc_matrix(3,1+2*n)   = 1; bc_target(3) = 0;         % I3(za) = 0

bc_matrix(4,n)       = 1; bc_target(4) = ia(1);     % I1(zb) = 0;
bc_matrix(5,2*n)     = 1; bc_target(5) = 0;         % I2(zb) = calculated;
bc_matrix(6,3*n)     = 1; bc_target(6) = ia(2);     % I3(zb) = 0;

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; gradient_constraint; bc_target];
constraint_matrix = sparse(constraint_matrix);

% Initial guess is the simple linear solution
i1_guess = linspace(i0(1),ia(1),n);
i2_guess = linspace(i0(2),0,n);
i3_guess = linspace(0,ia(2),n);
init_guess = [i1_guess i2_guess i3_guess];

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2);

 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,[],[],[],options);
 
i1_0a = x(n1);
i2_0a = x(n2);
i3_0a = x(n3);

G_0a = i1_0a.*g1 + i2_0a.*g2 + i3_0a.*g3;
B_0a = i1_0a.*b1 + i2_0a.*b2 + i3_0a.*g3;

%% Plot it

figure(2);
clf
plot(1e3*z_0a,i1_0a,'-','linewidth',1','color',co(1,:))
hold on
plot(1e3*z_0a,i2_0a,'-','linewidth',1','color',co(2,:))
plot(1e3*z_0a,i3_0a,'-','linewidth',1','color',co(3,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[z_init za]);
ylim([-50 70]);

%% Zone a to b calculation
n = 100;                % Points in this zone to evaluate
z_ab = linspace(za,zb,n);  % Vector of points to calculate

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);
n4 = (1+3*n):(4*n);

% Recalculate the field at this specific mesh because it is numerically
% easier to solve the problem on a mesh that is equally spaced between the
% high symmetry points.

% Find the Bfield at these points
b1 = interp1(Z,B(:,1),z_ab,'spline');
b2 = interp1(Z,B(:,2),z_ab,'spline');
b3 = interp1(Z,B(:,3),z_ab,'spline');
b4 = interp1(Z,B(:,4),z_ab,'spline');

% Find the gradient at these points
g1 = interp1(Z,G(:,1),z_ab,'spline');
g2 = interp1(Z,G(:,2),z_ab,'spline');
g3 = interp1(Z,G(:,3),z_ab,'spline');
g4 = interp1(Z,G(:,4),z_ab,'spline');

% Construct the field constraint matrix
field_constraint_matrix = [diag(b1) diag(b2) diag(b3) diag(b4)];
field_constraint = ones(n,1)*B0;

% Construct the gradient constraint matrix
gradient_constraint_matrix = [diag(g1) diag(g2) diag(g3) diag(g4)];
gradient_constraint = ones(n,1)*G0;

% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,n*4);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = ia(1);     % I1(za) = calculated
bc_matrix(2,1+n)     = 1; bc_target(2) = 0;         % I2(za) = 0
bc_matrix(3,1+2*n)   = 1; bc_target(3) = ia(2);     % I3(za) = calculated
bc_matrix(4,1+3*n)   = 1; bc_target(4) = 0;         % I4(za) = 0;

bc_matrix(5,n)       = 1; bc_target(5) = 0;         % I1(zb) = 0;
bc_matrix(6,2*n)     = 1; bc_target(6) = ib(1);     % I2(zb) = calculated;
bc_matrix(7,3*n)     = 1; bc_target(7) = 0;         % I3(zb) = 0;
bc_matrix(8,4*n)     = 1; bc_target(8) = ib(2);     % I4(zb) = calculated;

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; gradient_constraint; bc_target];
constraint_matrix = sparse(constraint_matrix);

% Initial guess is the simple linear solution
i1_guess = linspace(ia(1),0,n);
i2_guess = linspace(0,ib(1),n);
i3_guess = linspace(ia(2),0,n);
i4_guess = linspace(0,ib(2),n);
init_guess = [i1_guess i2_guess i3_guess i4_guess];

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2) + ...
    sum(diff(curr(n4)).^2);

 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,[],[],[],options);
 
i1_ab = x(n1);
i2_ab = x(n2);
i3_ab = x(n3);
i4_ab = x(n4);

G_ab = i1_ab.*g1 + i2_ab.*g2 + i3_ab.*g3 + i4_ab.*g4;
B_ab = i1_ab.*b1 + i2_ab.*b2 + i3_ab.*b3 + i4_ab.*b4;

%% Plot it

figure(3);
clf
plot(1e3*z_ab,i1_ab,'-','linewidth',1','color',co(1,:))
hold on
plot(1e3*z_ab,i2_ab,'-','linewidth',1','color',co(2,:))
plot(1e3*z_ab,i3_ab,'-','linewidth',1','color',co(3,:))
plot(1e3*z_ab,i4_ab,'-','linewidth',1','color',co(4,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[za zb]);
ylim([-50 70]);
% yyaxis right
% plot(1e3*z_bc,G_ab,'k-')
% ylabel('gradient (G/cm)');
%% Zone b to c calculation
n = 100;                % Points in this zone to evaluate
z_bc = linspace(zb,zc,n);  % Vector of points to calculate

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);
n4 = (1+3*n):(4*n);

% Recalculate the field at this specific mesh because it is numerically
% easier to solve the problem on a mesh that is equally spaced between the
% high symmetry points.

% Find the Bfield at these points
b2 = interp1(Z,B(:,2),z_bc,'spline');
b3 = interp1(Z,B(:,3),z_bc,'spline');
b4 = interp1(Z,B(:,4),z_bc,'spline');
b5 = interp1(Z,B(:,5),z_bc,'spline');

% Find the gradient at these points
g2 = interp1(Z,G(:,2),z_bc,'spline');
g3 = interp1(Z,G(:,3),z_bc,'spline');
g4 = interp1(Z,G(:,4),z_bc,'spline');
g5 = interp1(Z,G(:,5),z_bc,'spline');

% Construct the field constraint matrix
field_constraint_matrix = [diag(b2) diag(b3) diag(b4) diag(b5)];
field_constraint = ones(n,1)*B0;

% Construct the gradient constraint matrix
gradient_constraint_matrix = [diag(g2) diag(g3) diag(g4) diag(g5)];
gradient_constraint = ones(n,1)*G0;

% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,n*4);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = ib(1);     % I2(za) = calculated
bc_matrix(2,1+n)     = 1; bc_target(2) = 0;         % I3(za) = 0
bc_matrix(3,1+2*n)   = 1; bc_target(3) = ib(2);     % I4(za) = calculated
bc_matrix(4,1+3*n)   = 1; bc_target(4) = 0;         % I5(za) = 0;

bc_matrix(5,n)       = 1; bc_target(5) = 0;         % I2(zb) = 0;
bc_matrix(6,2*n)     = 1; bc_target(6) = ic(1);     % I4(zb) = calculated;
bc_matrix(7,3*n)     = 1; bc_target(7) = 0;         % I4(zb) = 0;
bc_matrix(8,4*n)     = 1; bc_target(8) = ic(2);     % I5(zb) = calculated;

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; gradient_constraint; bc_target];
constraint_matrix = sparse(constraint_matrix);

% Initial guess is the simple linear solution
i2_guess = linspace(ib(1),0,n);
i3_guess = linspace(0,ic(1),n);
i4_guess = linspace(ib(2),0,n);
i5_guess = linspace(0,ic(2),n);

init_guess = [i2_guess i3_guess i4_guess i5_guess];

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2) + ...
    sum(diff(curr(n4)).^2);

 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,[],[],[],options);
 
i2_bc = x(n1);
i3_bc = x(n2);
i4_bc = x(n3);
i5_bc = x(n4);

G_bc = i2_bc.*g2 + i3_bc.*g3 + i4_bc.*g4 + i5_bc.*g5;
B_bc = i2_bc.*b2 + i3_bc.*b3 + i4_bc.*b4 + i5_bc.*g5;

%% Plot it

figure(4);
clf
plot(1e3*z_bc,i2_bc,'-','linewidth',1','color',co(2,:))
hold on
plot(1e3*z_bc,i3_bc,'-','linewidth',1','color',co(3,:))
plot(1e3*z_bc,i4_bc,'-','linewidth',1','color',co(4,:))
plot(1e3*z_bc,i5_bc,'-','linewidth',1','color',co(5,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[zb zc]);
ylim([-50 70]);
% yyaxis right
% plot(1e3*z_bc,G_bc,'k-')
% ylabel('gradient (G/cm)');

%% Zone c to d calculation
n = 100;                % Points in this zone to evaluate
z_cd = linspace(zc,zd,n);  % Vector of points to calculate

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);
n4 = (1+3*n):(4*n);

% Recalculate the field at this specific mesh because it is numerically
% easier to solve the problem on a mesh that is equally spaced between the
% high symmetry points.

% Find the Bfield at these points
b3 = interp1(Z,B(:,3),z_cd,'spline');
b4 = interp1(Z,B(:,4),z_cd,'spline');
b5 = interp1(Z,B(:,5),z_cd,'spline');
b6 = interp1(Z,B(:,6),z_cd,'spline');

% Find the gradient at these points
g3 = interp1(Z,G(:,3),z_cd,'spline');
g4 = interp1(Z,G(:,4),z_cd,'spline');
g5 = interp1(Z,G(:,5),z_cd,'spline');
g6 = interp1(Z,G(:,6),z_cd,'spline');

% Construct the field constraint matrix
field_constraint_matrix = [diag(b3) diag(b4) diag(b5) diag(b6)];
field_constraint = ones(n,1)*B0;

% Construct the gradient constraint matrix
gradient_constraint_matrix = [diag(g3) diag(g4) diag(g5) diag(g6)];
gradient_constraint = ones(n,1)*G0;

% Construct the boundary condition constraint matrix
bc_matrix = zeros(8,n*4);
bc_target = zeros(8,1);
bc_matrix(1,1)       = 1; bc_target(1) = ic(1);     % I2(za) = calculated
bc_matrix(2,1+n)     = 1; bc_target(2) = 0;         % I3(za) = 0
bc_matrix(3,1+2*n)   = 1; bc_target(3) = ic(2);     % I4(za) = calculated
bc_matrix(4,1+3*n)   = 1; bc_target(4) = 0;         % I5(za) = 0;

bc_matrix(5,n)       = 1; bc_target(5) = 0;         % I2(zb) = 0;
bc_matrix(6,2*n)     = 1; bc_target(6) = id(1);     % I4(zb) = calculated;
bc_matrix(7,3*n)     = 1; bc_target(7) = 0;         % I4(zb) = 0;
bc_matrix(8,4*n)     = 1; bc_target(8) = id(2);     % I5(zb) = calculated;

% Assemble all constraints
constraint_matrix = [field_constraint_matrix; gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; gradient_constraint; bc_target];
constraint_matrix = sparse(constraint_matrix);

% Initial guess is the simple linear solution
i3_guess = linspace(ic(1),0,n);
i4_guess = linspace(0,id(1),n);
i5_guess = linspace(ic(2),0,n);
i6_guess = linspace(0,id(2),n);

init_guess = [i3_guess i4_guess i5_guess i6_guess];

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2) + ...
    sum(diff(curr(n4)).^2);

 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,[],[],[],options);
 
i3_cd = x(n1);
i4_cd = x(n2);
i5_cd = x(n3);
i6_cd = x(n4);

G_cd = i3_cd.*g3 + i4_cd.*g4 + i5_cd.*g5 + i6_cd.*g6;
B_cd = i3_cd.*b3 + i4_cd.*b4 + i5_cd.*b5 + i6_cd.*b6;

%% Plot it

figure(5);
clf
plot(1e3*z_cd,i3_cd,'-','linewidth',1','color',co(3,:))
hold on
plot(1e3*z_cd,i4_cd,'-','linewidth',1','color',co(4,:))
plot(1e3*z_cd,i5_cd,'-','linewidth',1','color',co(5,:))
plot(1e3*z_cd,i6_cd,'-','linewidth',1','color',co(6,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[zc zd]);
ylim([-50 70]);
% yyaxis right
% plot(1e3*z_cd,G_cd,'k-')
% ylabel('gradient (G/cm)');

%% Zone d to final calculation
n =100;                % Points in this zone to evaluate
z_dfinal = linspace(zd,z_final,n);  % Vector of points to calculate

n1 = (1):(n);
n2 = (1+n):(2*n);
n3 = (1+2*n):(3*n);

% Recalculate the field at this specific mesh because it is numerically
% easier to solve the problem on a mesh that is equally spaced between the
% high symmetry points.

% Find the Bfield at these points
b4 = interp1(Z,B(:,4),z_dfinal,'spline');
b5 = interp1(Z,B(:,5),z_dfinal,'spline');
b6 = interp1(Z,B(:,6),z_dfinal,'spline');

% Find the gradient at these points
g4 = interp1(Z,G(:,4),z_dfinal,'spline');
g5 = interp1(Z,G(:,5),z_dfinal,'spline');
g6 = interp1(Z,G(:,6),z_dfinal,'spline');

% Construct the field constraint matrix
field_constraint_matrix = [diag(b4) diag(b5) diag(b6)];
field_constraint = ones(n,1)*B0;

% Construct the gradient constraint matrix
gradient_constraint_matrix = [diag(g4) diag(g5) diag(g6)];
gradient_constraint = ones(n,1)*G0;


% Construct the boundary condition constraint matrix
bc_matrix = zeros(6,n*3);
bc_target = zeros(6,1);
bc_matrix(1,1)       = 1; bc_target(1) = id(1);     % I1(za) = calculated
bc_matrix(2,1+n)     = 1; bc_target(2) = 0;         % I2(za) = calculated
bc_matrix(3,1+2*n)   = 1; bc_target(3) = id(2);     % I3(za) = 0
bc_matrix(4,n)       = 1; bc_target(4) = 0;             % I1(zb) = 0;
bc_matrix(5,2*n)     = 1; bc_target(5) = ifinal(1);     % I2(zb) = calculated;
bc_matrix(6,3*n)     = 1; bc_target(6) = ifinal(2);     % I3(zb) = 0;
% Assemble all constraints
constraint_matrix = [field_constraint_matrix; gradient_constraint_matrix; bc_matrix];
constraint_vector = [field_constraint; gradient_constraint; bc_target];
constraint_matrix = sparse(constraint_matrix);

% Initial guess is the simple linear solution
 i4_guess = linspace(id(1),0,n);
 i5_guess = linspace(0,ifinal(1),n);
 i6_guess = linspace(id(1),ifinal(2),n);
init_guess = [i4_guess i5_guess i6_guess];

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2);

 x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector,[],[],[],options);
% 
i4_dfinal = x(n1);
i5_dfinal = x(n2);
i6_dfinal = x(n3);

G_dfinal = i4_dfinal.*g4 + i5_dfinal.*g5 + i6_dfinal.*g6;
B_dfinal = i4_dfinal.*b4 + i5_dfinal.*b5 + i6_dfinal.*b6;

%% Plot it

figure(6);
clf
plot(1e3*z_dfinal,i4_dfinal,'-','linewidth',1','color',co(4,:))
hold on
plot(1e3*z_dfinal,i5_dfinal,'-','linewidth',1','color',co(5,:))
plot(1e3*z_dfinal,i6_dfinal,'-','linewidth',1','color',co(6,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[zd z_final]);
ylim([-50 70]);

%% Final Output
data = struct;
data.z0a = z_0a;
data.z0a = z_0a;

% Raw
Z; %Position in meters 1xN
G; % Gradient in G/cm Nx6
B; % Field in G Nx6

% 0 to a
z_0a;
i1_0a;
i2_0a;
i3_0a;

% a to b
z_ab;
i1_ab;
i2_ab;
i3_ab;
i4_ab;

% b to c
z_bc;
i2_bc;
i3_bc;
i4_bc;
i5_bc;

% c to d
z_cd;
i3_cd;
i4_cd;
i5_cd;
i6_cd;

% d to fina
z_dfinal;
i4_dfinal;
i5_dfinal;
i6_dfinal;

% Combine them
i1 = [i1_0a i1_ab zeros(1,length(z_bc)) zeros(1,length(z_cd)) zeros(1,length(z_dfinal))];
i2 = [i2_0a i2_ab i2_bc zeros(1,length(z_cd)) zeros(1,length(z_dfinal))];
i3 = [i3_0a i3_ab i3_bc i3_cd zeros(1,length(z_dfinal))];
i4 = [zeros(1,length(z_0a)) i4_ab i4_bc i4_cd i4_dfinal];
i5 = [zeros(1,length(z_0a)) zeros(1,length(z_ab)) i5_bc i5_cd i5_dfinal];
i6 = [zeros(1,length(z_0a)) zeros(1,length(z_ab)) zeros(1,length(z_bc)) i6_cd i6_dfinal];
zz = [z_0a z_ab z_bc z_cd z_dfinal];

% if doSave 
    save('transport_calcs.mat','i1','i2','i3','i4','i5','i6','zz','coils')
% end
%%

figure(10);
clf
co=get(gca,'colororder');

plot([1 1]*za*1e3,[-50 60],'k:'); hold on
plot([1 1]*zb*1e3,[-50 60],'k:'); 
plot([1 1]*zc*1e3,[-50 60],'k:'); 
plot([1 1]*zd*1e3,[-50 60],'k:'); 

hold on

p1=plot(zz*1e3,i1,'-','color',co(1,:),'linewidth',1);
hold on
p2=plot(zz*1e3,i2,'-','color',co(2,:),'linewidth',1);
p3=plot(zz*1e3,i3,'-','color',co(3,:),'linewidth',1);
p4=plot(zz*1e3,i4,'-','color',co(4,:),'linewidth',1);
p5=plot(zz*1e3,i5,'-','color',co(5,:),'linewidth',1);
p6=plot(zz*1e3,i6,'-','color',co(6,:),'linewidth',1);

Gsolve = interp1(Z,G(:,1),zz).*i1+interp1(Z,G(:,2),zz).*i2 + ...
    interp1(Z,G(:,3),zz).*i3+interp1(Z,G(:,4),zz).*i4 + ...
    interp1(Z,G(:,5),zz).*i5+interp1(Z,G(:,6),zz).*i6;

Bsolve = interp1(Z,B(:,1),zz).*i1+interp1(Z,B(:,2),zz).*i2 + ...
    interp1(Z,B(:,3),zz).*i3+interp1(Z,B(:,4),zz).*i4 + ...
    interp1(Z,B(:,5),zz).*i5+interp1(Z,B(:,6),zz).*i6;


legend([p1,p2,p3,p4,p5],{'12a','12b','13','14','15','16'},'orientation','horizontal',...
    'location','southeast');

Gstr = ['G = ' num2str(mean(Gsolve)) '\pm' num2str(std(Gsolve),3) ' G/cm'];

Bstr = ['\DeltaB < ' num2str(range(Bsolve)) ' G'];

str = [Gstr newline Bstr];
% i1 = [i1_0a i1_ab]
xlim([0 max(zz)*1e3]);
ylim([-50 60]);
xlabel('position (mm)');
ylabel('current (A)');



text(.01,.99,str,'units','normalized','verticalalignment','top');

% yyaxis right
% set(gca,'YColor','k');
% ylabel('Gradient');
% plot(zz,Gsolve);




%% activate warning
warning on