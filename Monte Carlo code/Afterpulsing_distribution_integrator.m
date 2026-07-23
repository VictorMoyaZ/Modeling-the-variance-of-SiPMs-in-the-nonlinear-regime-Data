function [C] = Afterpulsing_distribution_integrator(a,Tau,Trec,delay,U,U0,Uch)

% VM: This function integrate the afterpulsing distribution detailed in:
% Rosado, Jaime, and Salvador Hidalgo. "Characterization and modeling of crosstalk and afterpulsing in Hamamatsu silicon photomultipliers." 
% Journal of Instrumentation 10.10 (2015): P10031. 

% VM: As input parameters it takes the diferent SiPM constant that the
%     model requires
%
%     As output parameter it calculate the normalization constant of the
%     distribution

n_bin=10^6;
t=linspace(delay,10*Tau+delay,n_bin)';
gain=1-exp(-t./Trec);
Trigger=(1-exp(-(gain.*U-U0)/Uch))./(1-exp(-(U-U0)/Uch));
Y=(t.^a).*exp(-t./Tau).*Trigger;

C=trapz(t,Y);

