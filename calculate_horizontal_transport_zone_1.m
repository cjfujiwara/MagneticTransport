%% Specify Boundary and Coils
x0 = 0;
x1 = 0.035;

alpha0 = 1;
alpha1 = 2.5;
G0 = 100;

i1 = 1;
i2 = 2;
i3 = 3;
n = 100;

xq = linspace(x0,x1,n);

%% Get Field

% Currents
I_mat_zone_1 = zeros(3,n);


% Field
B = zeros(1,n);

% Aspect Ratios
A  = linspace(alpha0,alpha1,n);

A = interp1([0  0.034 0.035],[1 2.5 2.5],xq,'linear');

A = interp1([0  0.01 0.03 0.035],[1 1.6 2.5 2.5],xq,'pchip');

% 
% x2alpha = @(x) (alpha0-alpha1)/x1^2*(x-x1).^2+alpha1;
% A = x2alpha(xq);

% Vertical Field Gradient
Gv  = ones(1,n)*G0;

% Horizontal Field Gradient
Gx  = -G0./(1+A);
Gy = zeros(1,n);

% Fields
B1 = interp1(X,Bx_all(i1,:),xq,'spline');
B2 = interp1(X,Bx_all(i2,:),xq,'spline');
B3 = interp1(X,Bx_all(i3,:),xq,'spline');

% Horizontal Gradient
Gx1 = interp1(X,Gx_all(i1,:),xq,'spline');
Gx2 = interp1(X,Gx_all(i2,:),xq,'spline');
Gx3 = interp1(X,Gx_all(i3,:),xq,'spline');

% Horizontal Gradient
Gy1 = interp1(X,Gy_all(i1,:),xq,'spline');
Gy2 = interp1(X,Gy_all(i2,:),xq,'spline');
Gy3 = interp1(X,Gy_all(i3,:),xq,'spline');

% Vertical Gradients
Gz1 = interp1(X,Gz_all(i1,:),xq,'spline');
Gz2 = interp1(X,Gz_all(i2,:),xq,'spline');
Gz3 = interp1(X,Gz_all(i3,:),xq,'spline');

for kk=1:length(xq)
    M = [B1(kk) B2(kk) B3(kk);
        Gz1(kk) Gz2(kk) Gz3(kk);
        Gx1(kk) Gx2(kk) Gx3(kk)];
    b = [B(kk);Gv(kk);Gx(kk)];
    I_mat_zone_1(:,kk) = mldivide(M,b);  
    Gy(kk) = Gy1(kk)*I_mat_zone_1(1,kk)+Gy2(kk)*I_mat_zone_1(2,kk)+Gy3(kk)*I_mat_zone_1(3,kk);
end

%% Plot It
hF = figure(11);
clf
hF.Color='w';
hF.Position=[10 60 350 850];
subplot(4,1,[1 2 3]);
p1=plot(xq*1e3,I_mat_zone_1(1,:),'color',co(i1,:),'linewidth',2);hold on;
p2=plot(xq*1e3,I_mat_zone_1(2,:),'color',co(i2,:),'linewidth',2);hold on;
p3=plot(xq*1e3,I_mat_zone_1(3,:),'color',co(i3,:),'linewidth',2);hold on;

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
