%% Properties

alpha = 1.707;
G0 = 100;
x_symm  = [0.068 0.1015 0.1310 0.1645 0.1940 0.2275 0.257 0.290];
n = 100;
ia = 3;
ib = 9;
jj = 1;

I_mat_zone_3 = zeros(length(coils),n*(ib-ia+1));
XQ = zeros(1,n*(ib-ia+1));
for ii = ia:ib
    i1 = ii;
    i2 = ii+1;
    i3 = ii+2;
    xq = linspace(x_symm(jj),x_symm(jj+1),n);

    nthis = ((jj-1)*n+1):jj*n;
    
    %Field
    B = zeros(1,n);    
    % Aspect Ratios
    A  = ones(1,n)*alpha; 

    % Vertical Field Gradient
    Gv  = ones(1,n)*G0;    
    % Horizontal Field Gradient
    Gx  = -G0./(1+A);
    Gy = zeros(1,n);

    % Currents
    imat = zeros(3,n);

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
        imat(:,kk) = mldivide(M,b);  
        Gy(kk) = Gy1(kk)*imat(1,kk)+Gy2(kk)*imat(2,kk)+Gy3(kk)*imat(3,kk);
    end

    I_mat_zone_3(i1,nthis)=imat(1,:);
    I_mat_zone_3(i2,nthis)=imat(2,:);
    I_mat_zone_3(i3,nthis)=imat(3,:);
    XQ(nthis) = xq;
   

    jj = jj+1;
end

%%  Plot It
hF = figure(13);
clf
hF.Color='w';
hF.Position = [10+700 60 350 850];
clear ps
ps=[];
subplot(4,1,[1 2 3]);
for ii = ia:(ib+2)
    ps(end+1)=plot(XQ*1e3,I_mat_zone_3(ii,:),'color',co(ii,:),'linewidth',2);hold on;
end
xlabel('position (mm)');
ylabel('current (A)');
set(gca,'fontsize',8,'box','on','linewidth',1);
ylim([0 100]);

xlim([XQ(1) XQ(end)]*1e3);


legend(ps,{coils(ia:(ib+2)).Name},'location','northwest');

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
xlim([xq(1) xq(end)]*1e3);
xlabel('position (mm)');
