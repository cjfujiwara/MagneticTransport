%% Specify Boundary and Coils

% Position
x0 = 0.325;                 % Starting Position
x1 = 0.360;                 % Ending Position
n = 200;                    % Number of position points
xq5 = linspace(x0,x1,n);    % Position Vector

% Aspect Ratios
alpha0 = 1.707;
alpha1 = 1;
A  = linspace(alpha0,alpha1,n);

% Aspect ratio interpolate
A = interp1(...
    [x0 0.335 0.345 0.350 0.355 0.360;],...
    [1.707 1.69 1.57 1.44 1.24 1],...
    xq5,'pchip');

% Target field gradient
G0 = 100;

% Coil indeces
i1 = 11;
i2 = 12;
i3 = 13;

%% Get Field

% Currents
I_mat_zone_5 = zeros(3,n);


% Field
B = zeros(1,n);

% Vertical Field Gradient
Gv  = ones(1,n)*G0;

% Horizontal Field Gradient
Gx  = -G0./(1+A);
Gy = zeros(1,n);

% Fields
B1 = interp1(X,Bx_all(i1,:),xq5,'spline');
B2 = interp1(X,Bx_all(i2,:),xq5,'spline');
B3 = interp1(X,Bx_all(i3,:),xq5,'spline');

% Horizontal Gradient
Gx1 = interp1(X,Gx_all(i1,:),xq5,'spline');
Gx2 = interp1(X,Gx_all(i2,:),xq5,'spline');
Gx3 = interp1(X,Gx_all(i3,:),xq5,'spline');

% Horizontal Gradient
Gy1 = interp1(X,Gy_all(i1,:),xq5,'spline');
Gy2 = interp1(X,Gy_all(i2,:),xq5,'spline');
Gy3 = interp1(X,Gy_all(i3,:),xq5,'spline');

% Vertical Gradients
Gz1 = interp1(X,Gz_all(i1,:),xq5,'spline');
Gz2 = interp1(X,Gz_all(i2,:),xq5,'spline');
Gz3 = interp1(X,Gz_all(i3,:),xq5,'spline');
Bout = zeros(1,n);

for kk=1:length(xq5)
    M = [B1(kk) B2(kk) B3(kk);
        Gz1(kk) Gz2(kk) Gz3(kk);
        Gx1(kk) Gx2(kk) Gx3(kk)];
    b = [B(kk);Gv(kk);Gx(kk)];
    I_mat_zone_5(:,kk) = mldivide(M,b);  
    Gy(kk) = Gy1(kk)*I_mat_zone_5(1,kk)+Gy2(kk)*I_mat_zone_5(2,kk)+Gy3(kk)*I_mat_zone_5(3,kk);

    Bout(kk) = B1(kk)*I_mat_zone_5(1,kk)+B2(kk)*I_mat_zone_5(2,kk)+B3(kk)*I_mat_zone_5(3,kk);
end

%% Plot It
hF = figure(15);
clf
hF.Color='w';
hF.Position=[1400+10 60 350 850];
subplot(4,1,[1 2 3]);
p1=plot(xq5*1e3,I_mat_zone_5(1,:),'color',co(i1,:),'linewidth',2);hold on;
p2=plot(xq5*1e3,I_mat_zone_5(2,:),'color',co(i2,:),'linewidth',2);hold on;
p3=plot(xq5*1e3,I_mat_zone_5(3,:),'color',co(i3,:),'linewidth',2);hold on;

xlabel('position (mm)');
ylabel('current (A)');
set(gca,'fontsize',8,'box','on','linewidth',1);
ylim([0 100]);

xlim([xq5(1) xq5(end)]*1e3)


legend([p1 p2 p3],legStr([i1 i2 i3]),'location','northwest');

subplot(4,1,4);
c = get(gca,'colororder');
yyaxis left
plot(xq5*1e3,Gv,'color',c(1,:),'linewidth',2);
ylim([90 110]);
ylabel('gradient (G/cm)');
xlabel('position (mm)');
yyaxis right
plot(xq5*1e3,Gy./Gx,'color',c(2,:),'linewidth',2);hold on;
set(gca,'fontsize',8,'box','on','linewidth',1);
ylabel('aspect ratio');
ylim([.9 3]);
xlim([xq5(1) xq5(end)]*1e3)
