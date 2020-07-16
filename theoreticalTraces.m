T=1;
xi=0;
xf=1;

tVec=linspace(0,T,100);


% Theoretical curve which minizmies the jerk from Jervis thesis
x=@(t) xi+(xf-xi)*(10*(t/T).^3-15*(t/T).^4+6*(t/T).^5);

s=.2*T;

figure(20);
set(gcf,'color','w');
clf
plot(tVec,x(tVec),'linewidth',2)
hold on
plot(tVec,xf*(1/2)*(1+erf((tVec-T/2)/(sqrt(2)*s))),'linewidth',2)


xlabel('time');
ylabel('position');

legend({'min jerk','erf'});


set(gca,'fontsize',14,'box','on','linewidth',1,'fontname','times');