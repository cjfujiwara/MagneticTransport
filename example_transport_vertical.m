%% Construct the coils

% Solve vertical transport in the "easy" regime
coils = makeVerticalCoils;

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
    z_centers(ii) = coils(ii).zbot+coils(2).Height/2;
end

% Idenfity High Symmetry Points
z_init = 0;
za = z_centers(2);
zb = z_centers(3);
zc = z_centers(4);
zd = z_centers(5);
z_final = z_centers(5)+z_centers(6);

z_symmetry = [z_init za zb zc zd z_final];
%% Initialize data vectors
N = 1e4; % Number of positions to evalulate
dL  = 1e-4; % Distance separation for calculating the gradient
Z = linspace(z_init,z_final,N);

B0 = 0;     % Target field in Gauss
G0 = 100;     % Target gradient in Gauss/cm

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
legend({'12a','12b','13','14','15','16'})

subplot(212)
for kk=1:length(coils)
    plot(Z*1e3,G(:,kk),'-','linewidth',1,'color',co(kk,:));
    hold on
end
xlabel('position (mm)')
ylabel('dB/dz @ 1A (G/cm)')
legend({'12a','12b','13','14','15','16'})


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

x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector);

i1_ab = x(n1);
i2_ab = x(n2);
i3_ab = x(n3);
i4_ab = x(n4);

Gab = i1_ab.*g1 + i2_ab.*g2 + i3_ab.*g3 + i4_ab.*g4;
Bab = i1_ab.*b1 + i2_ab.*b2 + i3_ab.*b3 + i4_ab.*b4;

%% Plot it

figure(2);
clf
plot(1e3*z_ab,i1_ab,'-','linewidth',1','color',co(1,:))
hold on
plot(1e3*z_ab,i2_ab,'-','linewidth',1','color',co(2,:))
plot(1e3*z_ab,i3_ab,'-','linewidth',1','color',co(3,:))
plot(1e3*z_ab,i4_ab,'-','linewidth',1','color',co(4,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[za zb]);
ylim([-50 50]);

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

x = fmincon(func,init_guess,[],[],constraint_matrix,constraint_vector);

i2_bc = x(n1);
i3_bc = x(n2);
i4_bc = x(n3);
i5_bc = x(n4);

Gbc = i2_bc.*g2 + i3_bc.*g3 + i4_bc.*g4 + i5_bc.*g5;
Bbc = i2_bc.*b2 + i3_bc.*b3 + i4_bc.*b4 + i5_bc.*g5;

%% Plot it

figure(3);
clf
plot(1e3*z_bc,i2_bc,'-','linewidth',1','color',co(2,:))
hold on
plot(1e3*z_bc,i3_bc,'-','linewidth',1','color',co(3,:))
plot(1e3*z_bc,i4_bc,'-','linewidth',1','color',co(4,:))
plot(1e3*z_bc,i5_bc,'-','linewidth',1','color',co(5,:))
xlabel('position (mm)');
ylabel('current (A)');
xlim(1e3*[zb zc]);
ylim([-50 50]);



















%% Zone 1 Solution OLD
zvec = linspace(za,zb,100);

c1 = coils(1).Coil;
c2 = coils(2).Coil;
c3 = coils(3).Coil;
c4 = coils(4).Coil;
dL  = 0.0001;

% Gradient and field value
G = 100;
B = 0;
b = [B; G];

% Field at BC1
[~,~,Bz1]=fieldCoil_3D(0,0,za,c1);      % Coil 1
[~,~,Bz2]=fieldCoil_3D(0,0,za,c2);      % Coil 2
[~,~,Bz3]=fieldCoil_3D(0,0,za,c3);      % Coil 3
[~,~,Bz4]=fieldCoil_3D(0,0,za,c4);      % Coil 4
Bz1 = Bz1*1e4;
Bz2 = Bz2*1e4;
Bz3 = Bz3*1e4;         
Bz4 = Bz4*1e4;

% Gradient at BC1
[~,~,Bz1p]=fieldCoil_3D(0,0,za+dL,c1);
[~,~,Bz2p]=fieldCoil_3D(0,0,za+dL,c2);
[~,~,Bz3p]=fieldCoil_3D(0,0,za+dL,c3);
[~,~,Bz4p]=fieldCoil_3D(0,0,za+dL,c4);
[~,~,Bz1n]=fieldCoil_3D(0,0,za-dL,c1);
[~,~,Bz2n]=fieldCoil_3D(0,0,za-dL,c2);
[~,~,Bz3n]=fieldCoil_3D(0,0,za-dL,c3);
[~,~,Bz4n]=fieldCoil_3D(0,0,za-dL,c4);
Gz1 = (Bz1p-Bz1n)/(2*dL);Gz1 = Gz1*1e2;
Gz2 = (Bz2p-Bz2n)/(2*dL);Gz2 = Gz2*1e2;
Gz3 = (Bz3p-Bz3n)/(2*dL);Gz3 = Gz3*1e2;  
Gz4 = (Bz4p-Bz4n)/(2*dL);Gz4 = Gz4*1e2;

M1 = [Bz1 Bz3; Gz1 Gz3];
i1 = M1\b;
i1 = [i1(1) 0  i1(2) 0];


% Field at BC1
[~,~,Bz1]=fieldCoil_3D(0,0,zb,c1);      % Coil 1
[~,~,Bz2]=fieldCoil_3D(0,0,zb,c2);      % Coil 2
[~,~,Bz3]=fieldCoil_3D(0,0,zb,c3);      % Coil 3
[~,~,Bz4]=fieldCoil_3D(0,0,zb,c4);      % Coil 4
Bz1 = Bz1*1e4;
Bz2 = Bz2*1e4;
Bz3 = Bz3*1e4;         
Bz4 = Bz4*1e4;

% Gradient at BC1
[~,~,Bz1p]=fieldCoil_3D(0,0,zb+dL,c1);
[~,~,Bz2p]=fieldCoil_3D(0,0,zb+dL,c2);
[~,~,Bz3p]=fieldCoil_3D(0,0,zb+dL,c3);
[~,~,Bz4p]=fieldCoil_3D(0,0,zb+dL,c4);
[~,~,Bz1n]=fieldCoil_3D(0,0,zb-dL,c1);
[~,~,Bz2n]=fieldCoil_3D(0,0,zb-dL,c2);
[~,~,Bz3n]=fieldCoil_3D(0,0,zb-dL,c3);
[~,~,Bz4n]=fieldCoil_3D(0,0,zb-dL,c4);
Gz1 = (Bz1p-Bz1n)/(2*dL);Gz1 = Gz1*1e2;
Gz2 = (Bz2p-Bz2n)/(2*dL);Gz2 = Gz2*1e2;
Gz3 = (Bz3p-Bz3n)/(2*dL);Gz3 = Gz3*1e2;  
Gz4 = (Bz4p-Bz4n)/(2*dL);Gz4 = Gz4*1e2;

M2 = [Bz2 Bz4; Gz2 Gz4];
i2 = M2\b;
i2 = [0 i2(1) 0  i2(2)];

%% Calculate all fields
Bz = zeros(length(zvec),4);
Gz = zeros(length(zvec),4);
curr_lin = zeros(length(zvec),4);

for jj=1:length(zvec)
    z0 = zvec(jj);

    % Field Values
    [~,~,Bz1]=fieldCoil_3D(0,0,z0,c1);
    [~,~,Bz2]=fieldCoil_3D(0,0,z0,c2);
    [~,~,Bz3]=fieldCoil_3D(0,0,z0,c3);
    [~,~,Bz4]=fieldCoil_3D(0,0,z0,c4);

Bz1 = Bz1*1e4;
Bz2= Bz2*1e4;
Bz3 = Bz3*1e4;
Bz4 = Bz4*1e4;

    % Z Gradient
    [~,~,Bz1p]=fieldCoil_3D(0,0,z0+dL,c1);
    [~,~,Bz2p]=fieldCoil_3D(0,0,z0+dL,c2);
    [~,~,Bz3p]=fieldCoil_3D(0,0,z0+dL,c3);
    [~,~,Bz4p]=fieldCoil_3D(0,0,z0+dL,c4);
    [~,~,Bz1n]=fieldCoil_3D(0,0,z0-dL,c1);
    [~,~,Bz2n]=fieldCoil_3D(0,0,z0-dL,c2);
    [~,~,Bz3n]=fieldCoil_3D(0,0,z0-dL,c3);
    [~,~,Bz4n]=fieldCoil_3D(0,0,z0-dL,c4);

    Gz1 = (Bz1p-Bz1n)/(2*dL);Gz1 = Gz1*1e2;
    Gz2 = (Bz2p-Bz2n)/(2*dL);Gz2 = Gz2*1e2;
    Gz3 = (Bz3p-Bz3n)/(2*dL);Gz3 = Gz3*1e2;
    Gz4 = (Bz4p-Bz4n)/(2*dL);Gz4 = Gz4*1e2;
    Bz(jj,:) = [Bz1 Bz2 Bz3 Bz4];
    Gz(jj,:) = [Gz1 Gz2 Gz3 Gz4];

    curr1 = i1(1) + (i2(1)-i1(1))/(zb-za)*(z0-za);
    curr2 = i1(2) + (i2(2)-i1(2))/(zb-za)*(z0-za);
    curr3 = i1(3) + (i2(3)-i1(3))/(zb-za)*(z0-za);
    curr4 = i1(4) + (i2(4)-i1(4))/(zb-za)*(z0-za);


    curr_lin(jj,:) = [curr1 curr2 curr3 curr4];
end


%%

n1 = 1:length(zvec);
n2 = (length(zvec)+1):(2*length(zvec));
n3 = (2*length(zvec)+1):(3*length(zvec));
n4 = (3*length(zvec)+1):(4*length(zvec));

Ball = Bz(:);
Gall = Gz(:);
curr_lin_all = curr_lin(:);

Gset = G*ones(length(zvec),1);
Bset = zeros(length(zvec),1);

% Magnetic Field Constraint
A1 = zeros(length(zvec),4*length(zvec));
B1 = zeros(length(zvec),1);
for kk=1:length(zvec)
    A1(kk,kk) = Bz(kk,1);
    A1(kk,kk+length(zvec)) = Bz(kk,2);
    A1(kk,kk+2*length(zvec)) = Bz(kk,3);
    A1(kk,kk+3*length(zvec)) = Bz(kk,4);
    B1(kk) = 0;
end

% Field Gradient Constraint
A2 = zeros(length(zvec),4*length(zvec));
B2 = zeros(length(zvec),1);
for kk=1:length(zvec)
    A2(kk,kk) = Gz(kk,1);
    A2(kk,kk+length(zvec)) = Gz(kk,2);
    A2(kk,kk+2*length(zvec)) = Gz(kk,3);
    A2(kk,kk+3*length(zvec)) = Gz(kk,4);
    B2(kk) = 100;
end

% Boundary Condition Constrain
Abc = zeros(8,4*length(zvec));
Bbc = zeros(8,1);

Abc(1,1) = 1;Bbc(1) = i1(1);                
Abc(2,1+length(zvec)) = 1;Bbc(2) = i1(2);
Abc(3,1+2*length(zvec)) = 1;Bbc(3) = i1(3);
Abc(4,1+3*length(zvec)) = 1;Bbc(4) = i1(4);
Abc(5,length(zvec)) = 1;Bbc(5) = i2(1);                
Abc(6,2*length(zvec)) = 1;Bbc(6) = i2(2);                
Abc(7,3*length(zvec)) = 1;Bbc(7) = i2(3);                
Abc(8,4*length(zvec)) = 1;Bbc(8) = i2(4);                


A = [A1;A2;Abc];
B = [B1;B2;Bbc];
x0 = curr_lin_all;

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2) + ...
    sum(diff(curr(n4)).^2);

x = fmincon(func,x0,[],[],A,B);

Gfind = x(n1).*Gz(:,1)+x(n2).*Gz(:,2)+x(n3).*Gz(:,3)+x(n4).*Gz(:,4);
Bfind = x(n1).*Bz(:,1)+x(n2).*Bz(:,2)+x(n3).*Bz(:,3)+x(n4).*Bz(:,4);
%%
figure(10);
clf

subplot(121);
plot(zvec,x(n1));
hold on
plot(zvec,x(n2));
plot(zvec,x(n3));
plot(zvec,x(n4));

subplot(122);
plot(zvec,Gfind); 
hold on
plot(zvec,Bfind);

