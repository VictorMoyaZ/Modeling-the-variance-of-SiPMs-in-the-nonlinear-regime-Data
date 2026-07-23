function [I] = instant_current(t,t0,I0,trec)

% VM: This function calculates the instant current created by an event

% VM: As input parameter this function takes:
%     t: The time at the current will be evaluated
%     t0: The time at the event happened, it t<t0 the current will be 0
%     I0: The current at the time the event happened
%     trec: The characteristic recovery time of the SiPM

% VM: As output parameter it calculate the instant current at the specified
%     time

if t<t0
    I=0;
else
    I=I0*exp(-(t-t0)/trec);
end
end

