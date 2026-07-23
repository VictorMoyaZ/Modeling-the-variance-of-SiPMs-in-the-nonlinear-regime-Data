function [Nap_t] = get_delayed_noise_time(events,P_list,T_list)

%This routine provides the time stamps for the afterpulse photons. 

% VM: used to be get_ap_time by interpolation

% VM: as input parameters it takes 
%     events: the number of events
%     P_list: the list of integrated probability until a certain time 
%     T_list: the list of times until the probability distribution have
%             been integrated

% VM: as output parameter it generates a vector of randomly selected times

     Nap_t=zeros(events,1);
    for n=1:events
        Nap_t(n,1)=interp1(P_list,T_list,rand(1,1));
 
    end  
    
    