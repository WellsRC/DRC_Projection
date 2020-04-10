function [beta,sigma,M,M2,gamma,theta,delta,P,MortR]=ParameterOutput(Amin,R0E)
%% PARAMETEROUTPUT returns the parameters based on the number of age classes
%specified in the model
% Input
% Amin - The minimum age of the different classes
% pd - probability of death in hospital
% Output
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
%% Paramter specification
sigma=1/4;
delta=1/(2*(7-4));
gamma=1; % onset to fever https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30566-3/fulltext?fbclid=IwAR1828JEbKgFGZLY9khGrg2Do4sDgKCXaURjEU3EC_CYWUKNl7hz3QXtI3k
theta=1-[0.8;0.8;0.4;0.2]; % PNAS Seyed
[M,M2,P,MortR]=DemoDRC(Amin);
beta=CalcR0NH(R0E,P,sigma,delta,M);
end

