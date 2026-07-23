function [flag_selector] = flag_selector(data,flag)

% VM: This function extract all the informations of events with a certain
%     flag
%
%     As inpunt parameters it takes the data struct created by
%     sipm_sim_array.m and a certain flag
%
%     As output parameters it generates a data struct named flag_selector
%     that contains al the information of the events of a certain flags
%

flag_selector = struct('count',0,'t',[], 'flag',[],'i',[],'j',[],'gain',[],'analysed_count',0);
flag_selector.count=length(data.flag==flag);
flag_selector.t=data.t(data.flag==flag);
flag_selector.i=data.i(data.flag==flag);
flag_selector.j=data.j(data.flag==flag);
flag_selector.gain=data.gain(data.flag==flag);
flag_selector.analysed_count=length(data.flag==flag);
end

