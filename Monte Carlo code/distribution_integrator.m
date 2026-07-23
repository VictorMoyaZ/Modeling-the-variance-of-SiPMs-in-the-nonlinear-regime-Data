function [P,T,C] = distribution_integrator(a,Tau,N,delay)

% VM: This function integrates a distribution: Y(t)=C*t^a*exp(t/Tau) and
%     gives back C and the integrated distribution between a minumum delay 
%     (before this time the distribution is supposed to be 0) and a vector
%     of times. 


% VM: as input parameters
%     a: is the exponent of the distribution
%     Tau: is the characteristic time of the distribution
%     N: is the number of pointsof the output vectors
%     delay: is the minimun time of the distribution 

% VM: as output parameter it generates:
%     P: is the vector of the integrated probability per time
%     T: is the vector of times until the probability vector have been
%        integrated
%     C: is the normalization constant of the distribution

n_bin=10^3;
X=linspace(delay,10*Tau+delay,n_bin)';
Y=(X.^a).*exp(-X./Tau);
C=trapz(X,Y);



T=logspace(log10(delay),log10(10*Tau),N)';
P=zeros(N,1);
P(N)=1;
p=0;
i=1;
for n=2:(N-1)
    while X(i,1)<=T(n,1)
        y=Y(i,1);
        y_max=Y(i+1,1);
        dx=X(i+1,1)-X(i,1);
        p=p+dx*(y_max+y)/2/C;
        
        i=i+1;
    end
  
    P(n,1)=p;
end

[P, ia, ic] = unique(P,'stable');
T=T(ia);
