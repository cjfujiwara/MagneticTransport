% Solve vertical transport in the "easy" regime
coils = makeVerticalCoils;

% Find center points of inner coils
z=[];
for kk=2:(length(coils)-1)
    z(end+1) = coils(kk).zbot+coils(2).Height/2;
end
za = coils(2).zbot+coils(2).Height/2;
zb = coils(3).zbot+coils(3).Height/2;
zc = coils(4).zbot+coils(4).Height/2;
zd = coils(5).zbot+coils(5).Height/2;

%% Zone 1 Solution
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

    % Bme = curr1*Bz1 + curr2*Bz2 + curr3*Bz3 + curr4*Bz4;
    % Gme = curr1*Gz1 + curr2*Gz2 + curr3*Gz3 + curr4*Gz4;

    % Bcurve(jj) = Bme;
    % Gcurve(jj) = Gme;
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

A1 = zeros(length(zvec),4*length(zvec));
B1 = zeros(length(zvec),1);
for kk=1:length(zvec)
    A1(kk,kk) = Bz(kk,1);
    A1(kk,kk+length(zvec)) = Bz(kk,2);
    A1(kk,kk+2*length(zvec)) = Bz(kk,3);
    A1(kk,kk+3*length(zvec)) = Bz(kk,4);
    B1(kk) = 0;
end

A2 = zeros(length(zvec),4*length(zvec));
B2 = zeros(length(zvec),1);
for kk=1:length(zvec)
    A2(kk,kk) = Gz(kk,1);
    A2(kk,kk+length(zvec)) = Gz(kk,2);
    A2(kk,kk+2*length(zvec)) = Gz(kk,3);
    A2(kk,kk+3*length(zvec)) = Gz(kk,4);
    B2(kk) = 100;
end

% bbc = zeros(
% Abc = zeros(8,4*length(zvec));
% Abc(1,1) = 1;

% b = zeros(1,4*length(zvec));
% b(1) = i1(1);



A = [A1;A2];
B = [B1; B2];
x0 = curr_lin_all;

func = @(curr) sum(diff(curr(n1)).^2) + ...
    sum(diff(curr(n2)).^2) + ...
    sum(diff(curr(n3)).^2) + ...
    sum(diff(curr(n4)).^2);

x = fmincon(func,x0,[],[],A,B);
% diff(curr(1:n1)).^2+diff(curr(1:n1)).^2

