function [T_lim,Mean_gain,Events] = gain_evolution(gain,T_gain,T_max,Nbin,Zeros,T_medida_0,I_mean)

% VM: This function calculate the mean gain, number of events for specific
%     time bins, also it plots the mean gain, number of events and the
%     output current

% VM: as input parameters this function takes:
%     gain: vector of gains
%     T_gain: time asociated with the gain vector
%     T_max: maximum value of time bin
%     Nbin: number of time bins
%     Zeros: if equal to 1 unsuccessful events will be included in
%     calculation, otherwise they will be ignored
%     T_medida_0: The time vector asociated with the output current
%     I_mean: Output current vector

% VM: as output parameters:
%     T_lim: vector containing the time bins
%     Mean_gain: mean gain vector
%     Events: vector with the number of events

% VM: Generating the time bin limits
T_lim=linspace(0,T_max,Nbin+1);

% VM: if Zeros = 1 unsuccessful events will be included in the mean gain and event calculation 

% if Zeros==0
%     Gain=gain(gain>0);
%     T_Gain=T_gain(gain>0);
% else
%     Gain=gain;
%     T_Gain=T_gain;
% end

% VM: preparing data to be treated
Gain=gain;
T_Gain=T_gain;

[T_Gain, order] = sort(T_Gain,'ascend');
Gain=Gain(order);
N=length(Gain);
j=1;
A=[];
B=0;
Mean_gain=[];   
Events=[];

% VM: Obtaining the number of events and the mean gain in the time bin

for i=1:N
    if T_lim(j)<T_Gain(i)<T_lim(j+1)
        A=[A,Gain(i)];
        B=B+1;
    else
        if ~isempty(A)
            Mean_gain=[Mean_gain,mean(A)];
        else
            Mean_gain=[Mean_gain,0];
        end
        Events=[Events,B];
        A=[];
        B=0;
        j=j+1;
    end
    if T_lim(j) == T_max
        
        break
    end
end

% VM: Graphical representation of the mean gain, the number of events and the
%     output current

T_lim=T_lim(1:j-1);
subplot(3, 1, 1)
bar( T_lim*10^9, Mean_gain);
xlabel('Time ns')
ylabel('Mean Gain')
subplot(3, 1, 2)
bar( T_lim*10^9, Events)
xlabel('Time ns')
ylabel('Events')
subplot(3, 1, 3)
plot(T_medida_0*10^9, I_mean)
xlabel('Time ns')
ylabel('Mean current (A)')
end

