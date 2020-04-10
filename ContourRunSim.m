clear;
close all;
Amin=[0 20 50 65];
TS=100000;
A=length(Amin);
R0v=linspace(1.1,2.9,101);
Cases=zeros(4,101);
DeathsA=zeros(4,101);
DeathsB=zeros(4,101);
DeathsC=zeros(4,101);
DeathsD=zeros(4,101);
incf=linspace(0,0.05,51);
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
for mm=1:51
    DeathsA(mm,nn)=sum(YM1(end,CD).*(MortR+[incf(mm) 0 0 0]))./sum(P);
    DeathsB(mm,nn)=sum(YM1(end,CD).*(MortR +[0 incf(mm) 0 0]))./sum(P);
    DeathsC(mm,nn)=sum(YM1(end,CD).*(MortR +[0 0 incf(mm) 0]))./sum(P);
    DeathsD(mm,nn)=sum(YM1(end,CD).*(MortR+[0 0 0 incf(mm)]))./sum(P);
end
end
mm=max([max(DeathsA(:)) max(DeathsB(:)) max(DeathsC(:)) max(DeathsD(:))]);
mm=round(1000*mm)/1000;
load('Cmap.mat');
figure('units','normalized','outerposition',[0 0 1 1]);

contourf(R0v,100.*(MortR(1)+incf),DeathsA,'LineStyle','none');
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
ch=colorbar;
set(get(ch,'label'),'string','Proportion of population that will die from COVID-19','Rotation',270,'Position',[7.233650684356689,0.014676927892061,0]);
yh=ylabel('Case fatality ratio (%) for 0-19 age class','Fontsize',28);
box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.1:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
colormap(cmap);
caxis([0 mm]);
text(yh.Extent(1),max(ylim),'A','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure2A.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0 0 1 1]);

contourf(R0v,100.*(MortR(2)+incf),DeathsB,'LineStyle','none');
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
ch=colorbar;
set(get(ch,'label'),'string','Proportion of population that will die from COVID-19','Rotation',270,'Position',[7.233650684356689,0.014676927892061,0]);
yh=ylabel('Case fatality ratio (%) for 20-49 age class','Fontsize',28);
box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.1:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
colormap(cmap);
caxis([0 mm]);
text(yh.Extent(1),max(ylim),'B','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure2B.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0 0 1 1]);

contourf(R0v,100.*(MortR(3)+incf),DeathsC,'LineStyle','none');
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
ch=colorbar;
set(get(ch,'label'),'string','Proportion of population that will die from COVID-19','Rotation',270,'Position',[7.233650684356689,0.014676927892061,0]);
yh=ylabel('Case fatality ratio (%) for 50-64 age class','Fontsize',28);
box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.1:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
colormap(cmap);
caxis([0 mm]);
text(yh.Extent(1),max(ylim),'C','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure2C.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0 0 1 1]);

contourf(R0v,100.*(MortR(4)+incf),DeathsD,'LineStyle','none');
xlabel('Basic reproductive number ({\it R}_0)','Fontsize',28);
ch=colorbar;
set(get(ch,'label'),'string','Proportion of population that will die from COVID-19','Rotation',270,'Position',[7.233650684356689,0.014676927892061,0]);
yh=ylabel('Case fatality ratio (%) for 65+ age class','Fontsize',28);
box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[1.1:0.1:2.9],'Xminortick','on','Fontsize',26,'YminorTick','on');
colormap(cmap);
caxis([0 mm]);
text(yh.Extent(1),max(ylim),'D','Fontsize',38,'Fontweight','bold');
print(gcf,['Figure2D.png'],'-dpng','-r300');