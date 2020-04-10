% https://github.com/midas-network/COVID-19/tree/master/parameter_estimates/2019_novel_coronavirus
% Accessed April 9, 2020 (1:50 pm)
R0=[0.3
2.2
2.6
2.24
1.48
2.1
3.2
7.05
1.44
2.2
4.8
4.6
5.14
6.09
7.09
2.6
3.25
3.13
3.42
3.96
4.58
2
1.27
1.59
3.6
3.58
2.5
2.37
2.51
2.36
2.72
2.42
2.6
2.65
2.59
2.23
2.01
2.55
2.14
2.7
2.48
2.25
2.3
2.67
2.62
2.28
2.41
2.27
2.55
2.67
2.29
2.45
2.41
2.54
2.53
2.3
2.68
2.07
2.2
3.61
2.72
2.53
2.27
2.27
2.15
2.41
2.43
3.55
2.86
2.59
2.78
2.67
3.22
2.28
2.33
2.31
1.68
2.23
2.2
2.73
2.44
2.39
2.22
2.36
2.28
2.41
2.61
2.62
2.67
2.34
2.21
2.35
2.64
2.29
2.35
2.51
2.66
2.46
2.65
2.52
2.52
2.14
2.47
2.34
2.25
3.15
2.27
2.76
2.7
2.66
2.28
2.21
2.64
2.32
2.69
2.38
3.11
2.69
2.21
2.82
2.64
2.42
2.31
5.25
3.15
4.1
3.26
4.24
3.1
2.5
2.3
3.1
3.8
3
2.6
];
R0v=zeros(3,1);
R0v(1)=round(mean(R0),2);

SEM = std(R0)/sqrt(length(R0));               % Standard Error
ts = tinv([0.025  0.975],length(R0)-1);      % T-Score
R0v(2:3) = round(mean(R0) + ts*SEM,2);                      % Confidence Intervals

R0v
Amin=[0 20 50 65];
TS=365;
A=length(Amin);
Cases=zeros(1,3);
SCases=zeros(4,3);
Deaths=zeros(1,3);
SDeath=zeros(4,3);
for nn=1:3
R0E=R0v(nn);
[beta,sigma,M,M2,gamma,theta,delta,P,MortR]=ParameterOutput(Amin,R0E);
NState=8;
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
CS=7*A+[1:A]; % cumulative removal of infections

Cases(nn)=sum(YM1(end,CI));
SCases(:,nn)=(YM1(end,CI));
Deaths(nn)=sum(YM1(end,CD).*MortR);
SDeath(:,nn)=(YM1(end,CD).*MortR);
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
bar([1:4],SCases(:,1),'LineStyle','none','Facecolor',hex2rgb('#336B87'));
ylim([0 5*10^7]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:4],'XTicklabel',{'0-19','20-49','50-64','65+'},'Fontsize',18)
ylabel('Total number of infections','Fontsize',20);
xlabel('Age class','Fontsize',20);
box off;
xlim([0.5 4.5]);
subplot(2,2,2);
bar([1:4],SDeath(:,1),'LineStyle','none','Facecolor',hex2rgb('#336B87'));
ylim([0 2*10^5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:4],'XTicklabel',{'0-19','20-49','50-64','65+'},'Fontsize',18)
ylabel('Total number of deaths','Fontsize',20);
xlabel('Age class','Fontsize',20);
box off;
xlim([0.5 4.5]);

fprintf('Cases \n')
round(Cases)
fprintf('Deaths \n')
round(Deaths)

fprintf('Cases in 0-19 \n')
round(100.*SCases(1,:)./sum(SCases,1),2)

fprintf('Deaths in 50+ \n')
round(100.*(SDeath(3,:)+SDeath(4,:))./sum(SDeath,1),2)
