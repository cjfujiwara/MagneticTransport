c1 = dlmread('rev45coil1.txt' ,',',0,1);
c2 = dlmread('rev45coil2.txt' ,',',0,1);
c3 = dlmread('rev45coil3.txt' ,',',0,1);
c4 = dlmread('rev45coil4.txt' ,',',0,1);
c5 = dlmread('rev45coil5.txt' ,',',0,1);
c6 = dlmread('rev45coil6.txt' ,',',0,1);
%%

x = 1:1:174;

figure(4)
clf
co=get(gca,'colororder');
subplot(4,1,[1 2])
plot(x,c1,'-','color',co(1,:),'linewidth',1)
hold on
plot(x,c2,'-','color',co(2,:),'linewidth',1)
plot(x,c3,'-','color',co(3,:),'linewidth',1)
plot(x,c4,'-','color',co(4,:),'linewidth',1)
plot(x,c5,'-','color',co(5,:),'linewidth',1)
plot(x,c6,'-','color',co(6,:),'linewidth',1)

imat = [c1 c2 c3 c4 c5 c6];

xlabel('desired position (mm)')
ylabel('current (A)')
legend({'12a','12b','13','14','15','16'},'orientation','horizontal');
xlim([1 174])
ylim([-60 60]);
% B = B*1e4;
% G = G*1e2;

subplot(4,1,[3])
plot(x,g0vec/100,'k-')
ylim([85 110]);
xlabel('desired position (mm)')
ylabel('gradient (G/cm)')

subplot(4,1,[4])
plot(x,z0vec*1e3-x','k-')
 ylim([-6 2]);
xlabel('desired position (mm)')
ylabel('\Delta z (mm)')
%% Calculate Field and Gradient at symmetry points

z0vec = zeros(length(x),1);
g0vec = zeros(length(x),1);

for jj=1:length(x)
    b1 = B(:,1)*c1(jj);
    b2 = B(:,2)*c2(jj);
    b3 = B(:,3)*c3(jj);
    b4 = B(:,4)*c4(jj);
    b5 = B(:,5)*c5(jj);
    b6 = B(:,6)*c6(jj);

    Ball = b1+b2+b3+b4+b5+b6;
    [~,ind]=min(abs(Ball));
    z0vec(jj) = Z(ind);

    g0vec(jj) = (Ball(ind+1)-Ball(ind-1))/(2*(Z(2)-Z(1)));


end

figure(5)
subplot(121);
plot(x,z0vec,'-')
xlabel('desired positiion (mm)')
ylabel('calculate position (mm)')

subplot(122);
plot(x,g0vec,'-')
xlabel('desired positiion (mm)')
ylabel('calculate position (mm)')