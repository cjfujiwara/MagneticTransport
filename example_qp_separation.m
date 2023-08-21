% Solve vertical transport in the "easy" regime
warning off;

%% Calculation setting
N = 1e4; % Number of positions to evalulate
dL  = 1e-4; % Distance separation for calculating the gradient
Z = linspace(.170,.190,N);

B0 = 0;     % Target field in Gauss

% Vectors for magnetic fields and gradients
B = zeros(N,6);
G = zeros(N,6);

%% Initialize data vectors
coils = makeVerticalCoils;

coils = coils(5:6);

dz = linspace(0,25,100)*1e-3;
mt_center = zeros(length(dz),1);
mt_gradient = zeros(length(dz),1);
sep = zeros(length(dz),1);
for nn=1:length(dz)
    tic
    dz0=dz(nn);
    coils(end).zbot = 191.5e-3 + dz0; 
    coils(end).Coil = makeVerticalCoil(coils(end));
    sep(nn) = coils(end).zbot - coils(1).zbot;
    
    

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
    
    Bqp = B(:,2)-B(:,1);
    Gqp = G(:,2)-G(:,1);
    
    ind = find(diff(sign(Bqp)),1);
    
    mt_center(nn) = Z(ind);
    mt_gradient(nn) = Gqp(ind);
  toc
end

%%
ind = find(mt_gradient<7.39,1);
z_find = dz(ind);
sep_find = sep(ind);

center_find = mt_center(ind);


figure(200)
clf
co=get(gca,'colororder');
clf
cla
% 
% ax=axes;
yyaxis left
p1=plot(sep*1e3,mt_center*1e3,'-','linewidth',1,'color',co(1,:));hold on;
xlabel('separation (mm)');
ylabel('mt center (mm)');

plot(get(gca,'XLim'),[1 1]*center_find*1e3,'--','color',co(1,:));
set(gca,'YColor',co(1,:));

yyaxis right
p2 = plot(sep*1e3,mt_gradient,'k-','linewidth',1);
set(gca,'YColor','k');
ylabel('field gradient at B=0 (G/cm/A)');

plot([1 1]*sep_find*1e3,get(gca,'YLim'),'k--');

% ax.YAxis(2).Color = co(1,:);
str = ['\Deltaz = ' num2str(sep_find*1e3,3) ' mm (' num2str((191.5 - 133.7)) ' Yee)' ];
text(.05,.05,str,'units','normalized','fontsize',14,'verticalalignment','bottom');

legend([p1 p2],{'MT center','gradient'},'location','best');
% plot(

