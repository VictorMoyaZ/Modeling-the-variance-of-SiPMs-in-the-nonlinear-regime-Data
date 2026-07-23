function [Nph_t] = get_timestamp_photons(Nph, t,u,Poisson,Theta)

%This is a simple routine to evaluate the timestamp of the optical photons.
% The timestamps can also be obtained more rigorously using advanced
% software such as GEANT4. 

% VM: As input parameter this function takes:
%     Nph: Is the number of photons
%     t: Is the time asociated with each element of the integrated
%     probability vector
%     u: Is the integrated probability vector
%     Poisson: Is a flag, if it is equal to 1 it will randomly generate the
%     number of photons using a Poisson distribution with mean Nph

% VM: As output parameter it generates a vector containig the arrival time
%     of the incident photons

% VM: If Poisson equal to 1 it will randomly generate the number of
%     incident photons

    if Poisson == 1
        Nph = poissrnd(Nph,1); % Number of photons generated.
    else
        % Super-poisson case: Negative binomial
        mu=Nph;
        
        %var=mu/Theta;
        r=Theta*mu/(1-Theta);
        Nph = nbinrnd(r,Theta,1);
    end

% VM: It genrates Nph random numbers ant interpolate each corresponding
%     time in the t, u vectors 


    
    
    
    
    Nph_t = interp1(u,t,rand(Nph,1));
    
    
    