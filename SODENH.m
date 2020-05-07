function dxdt = SODENH(t,x,beta,sigma,M,M2,gamma,delta,theta,P,A)
%% dxdt = SODE()
% System of ODE to run the model for vaccination and contact tracing
%% Input parameters
%t - time variable
%x - State variable 
%     S=[1:A]; % Susceptible
%     E=A+[1:A]; % Incubation
%     I=2*A+[1:A]; %Symptomatic and will be hospitalized
%beta - probability of infection
%sigma - Rate from infection to symptoms
%tau - Contact tracing rate
%M - contact matrix (Size: AxA)
%M2 - contact matrix home (Size: AxA)
%gamma - rate to recovery symptomatic individual
%q - rate percentage of unvaccinated syptomatic case self-quaratine
%h - rate percentage of unvaccinated symptatic case being hospitalized
%delta - hospitalization rate
%mh - rate percentage for mortality in hospital
%mueH - mortality rate in hospital
%psiH - recover rate from hosptial
%mc - rate percentage for mortality in ICU
%mueC - mortality rate in ICU
%psiC - recover rate from ICU
%P -Population size
%A - Number of ages classes considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the vector and specify the index equations
dxdt=zeros(length(x),1);
% Index for dxdt and x to make readability of code easier
S=[1:A]; % Susceptible
E=A+[1:A]; % Incubation and non-vaccinated
I=2*A+[1:A]; % Incubation andvaccinated after infection
C=3*A+[1:A]; % Incubation andvaccinated after infection
Q=4*A+[1:A]; % Incubation andvaccinated after infection
CI=5*A+[1:A]; % cumulative infections
CD=6*A+[1:A]; % cumulative removal of infections

%% Susceptible and vaccinated population
dxdt(S)= -beta.*(M*((x(I)./P+x(C)./P))).*x(S)-beta.*(M2*(x(Q)./P)).*x(S);
%% Incubation period population
dxdt(E)= beta.*(M*(x(I)./P+x(C)./P)).*x(S)+beta.*(M2*(x(Q)./P)).*x(S)-sigma.*x(E);
%% Milder infection
dxdt(I)=(1-theta).*sigma.*x(E)-delta.*x(I);
%% More severe infection
dxdt(C)=theta.*sigma.*x(E)-gamma.*x(C);
%% At home 
dxdt(Q)=gamma.*x(C)-(delta*gamma)./(gamma-delta).*x(Q);
%% New infectious
dxdt(CI)=sigma.*x(E);

%% Resolved and infectious
dxdt(CD)=delta.*x(I)+(delta*gamma)./(gamma-delta).*x(Q);
end

