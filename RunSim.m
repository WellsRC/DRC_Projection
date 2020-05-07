clear;
close all;
Amin=[0 20 50 65];
TS=100000;
A=length(Amin);
R0v=linspace(1.2,2.9,101);
Cases=zeros(4,101);
Deaths=zeros(4,101);
for nn=1:101
R0E=R0v(nn);
[beta,sigma,M,M2,gamma,theta,delta,P,MortR]=ParameterOutput(Amin,R0E);
NState=7;
IC=zeros(NState*A,1);
IC(1:A)=P;
IC(1:A)=IC(1:A)-1;
IC(A+[1:A])=1;
options = odeset('RelTol',10^(-9),'AbsTol',(10^(-9).*ones(size(IC))),'NonNegative',1:(NState*A));
[TM1,YM1]=ode15s(@(t,y)SODENH(t,y,beta,sigma,M,M2,gamma,delta,theta,P,A),[0 TS],IC,options);

S=[1:A]; % Susceptible
E=A+[1:A]; % Incubation and non-vaccinated
I=2*A+[1:A]; % Incubation andvaccinated after infection
C=3*A+[1:A]; % Incubation andvaccinated after infection
Q=4*A+[1:A]; % Incubation andvaccinated after infection
CI=5*A+[1:A]; % cumulative infections
CD=6*A+[1:A]; % cumulative removal of infections

Cases(:,nn)=YM1(end,CI);
Deaths(:,nn)=YM1(end,CD).*MortR;
end

close all;

Colr=[hex2rgb('#E1B16A');
    hex2rgb('#78A5A3');    
    hex2rgb('#CE5A57');
    hex2rgb('#444C5C');];
figure('units','normalized','outerposition',[0 0 1 1]);
bb=area(R0v,Cases'./1000,'LineStyle','none');
for ii=1:4
 bb(ii).FaceColor=Colr(ii,:);
end
box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.2:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
xlim([1.2 2.9]);
ax = gca;
ax.YRuler.Exponent = 0;
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
legend({'Ages 0-19','Ages 20-49','Ages 50-64','Ages 65+'},'location','NorthWest','Fontsize',26);
legend boxoff;
yh=ylabel('Total number of cases (\times 1,000) ','Fontsize',28);
text(yh.Extent(1),max(ylim),'A','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure1A.png'],'-dpng','-r600');
saveas(gcf,'Test1','epsc')
figure('units','normalized','outerposition',[0 0 1 1]);
bb=area(R0v,Deaths'./1000,'LineStyle','none');
for ii=1:4
 bb(ii).FaceColor=Colr(ii,:);
end
box off;
xlim([1.2 2.9]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.2:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
ax = gca;
ax.YRuler.Exponent = 0;
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
yh=ylabel('Total number of deaths (\times 1,000)','Fontsize',28);
ylim([0 150]);
legend({'Ages 0-19','Ages 20-49','Ages 50-64','Ages 65+'},'location','NorthWest','Fontsize',26);
legend boxoff;
text(yh.Extent(1),max(ylim),'B','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure1B.png'],'-dpng','-r600');
saveas(gcf,'Test2','epsc')
