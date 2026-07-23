function [meas_sipm_par,fixed_sipm_par,crystal_par,flags,Nph,info] = read_input_array(sipm_no,V,pde,td,tr,ts,tt1,tt2,tc,Prob,Poisson,E)
%Read all the input parameters of all the SiPMs connected to the
%scintillator. Typically, to a single scintillator, multiple SiPMs are
%connected, so in this code, we read the measured parameters for all the SiPMs.
%A lot of symbols are the same as used and defined in the following paper:
%Jha, A. K., van Dam, H. T., Kupinski, M. A., & Clarkson, E. (2013).
%Simulating Silicon Photomultiplier Response to Scintillation Light.
%IEEE Transactions on Nuclear Science, 30(1), 336?351. http://doi.org/10.1109/TNS.2012.2234135


%Simulation flags

%VM: delayed crosstalk, trigger , simplified_readout flags added

flag_default_variation_parameters_overvoltage = 1;% This flag can be set to zero if
% the user wants to input specific variations with respect to overvoltage. 
% For example, if the user has measurements for how the different
% SiPM parameters, such as the afterpulsing rate, dark current etc. vary
% with the overvoltage and the user would like to use these measurements.

flags = struct('oct',1,'ap',0,'doct',0,'Biexponential_ap_distibution',0,'deadtime',0,'overvoltage',1,'trigger',1,'output_voltage',1,'dc',0,'partial_recovery',1,...
    'overvoltage_pde',0,'overvoltage_pde_meas',0,'overvoltage_pde_lin',0,...
    'overvoltage_dc',0,'overvoltage_dc_meas',0,'overvoltage_dc_lin',0,...
    'overvoltage_oct',0,'overvoltage_oct_meas',0,'overvoltage_oct_quad',0,...
    'overvoltage_ap',0,'overvoltage_ap_meas',0,'overvoltage_ap_quad',0,'simplified_readout',1);
% VM:

% oct:
% ap:
% doct:
% deadtime:
% overvoltage:
% trigger:
% output_voltage:
% dc:
% partial_recovery:
% overvoltage_pde:
% overvoltage_pde_meas:
% overvoltage_pde_lin:
% overvoltage_oct:
% overvoltage_oct_meas:
% overvoltage_oct_quad:
% overvoltage_ap: 
% overvoltage_ap_meas:
% overvoltage_ap_quad:
% simplified_readout:

% The parameters that describe the scintillation pulse


crystal_par = struct('Qsc',0,'Td',0,'Tr',0,'Ts',0,'Tt1',0,'Tt2',0,'Tc',0,'Poisson',0,'t',0,'u',0);
% VM:
% Qsc: scintillator yield
% Td: down time of exponential like light pulse
% Tr: rise time of exponential like light pulse
% Ts: characteristic time of a square light light pulse
% Tt1: characteristic rising time of a triangular light pulse
% Tt2: characteristic down time of a triangular light pulse
% Tc: characteristic time of a semicircular light pulse
% Poisson: flag of poisson fluctuations, if equal to 1 poisson fluctuations
%          will be activated
% t: time vector of the integrated probability vector. Used to randomly
%    select photon arrival times
% u: vector of integrated probability of the time distribution oh the photon arrival. Used to randomly
%    select photon arrival times 





%Td = 16e-9; % Decay time of the scintillator.
%Td = 1.24*10e-9;
Td=td;
%Td=0;
%Td=42*10e-9;

%Tr = 0.93e-9; % Rise time of the scintillator
Tr=tr;
%Tr = td/2;

%Ts = td;% VM: Time for square photon pulses
Ts=ts;

Tt1=tt1;
Tt2=tt2;

Tc=tc; % VM: Time for semicircular photon pulses

%Qsc = 70; % Photon yield for LaBR3
Qsc = 1;
%Qsc = 29;

crystal_par.Poisson = Poisson;
crystal_par.Qsc = Qsc;

crystal_par.Td = Td; %Decay time of scintillation pulse
crystal_par.Tr = Tr; %Rise time of scintillation pulse
crystal_par.Ts = Ts;
crystal_par.Tt1 = Tt1;
crystal_par.Tt2 = Tt2;
crystal_par.Tc = Tc;

%T_rec=10e-9; %VM: trec
%T_rec=29e-9;
T_rec=30e-9;
%T_rec=20e-9;

tmeas = 100*(Td+T_rec+Tr+Ts+Tt1+Tt2+Tc);
Nbar_ph = Qsc*E;
if Poisson == 0
        Nph =round(Nbar_ph);
        Nph = unique(Nph);
        Nph=Nph(Nph>0);
else 
    Nph =Nbar_ph;
end

if Tr==0 && Ts==0 && Td>0
        t = linspace(0,tmeas,1000);
        u = round(1 - exp(-t/Td),3);
    elseif Ts>0
        t = linspace(0,Ts,1000);   
        u=t/Ts;
    elseif Tr>0 && Td>0
        t = linspace(0,tmeas,1000);
        u = round(1 - Td/(Td - Tr)*exp(-t/Td) + Tr/(Td - Tr)* exp(-t/Tr),3);
    elseif Tt1>0 && Tt2==0
        t = linspace(0,Tt1,1000);
        u = (t/Tt1).^2 ; %rising triangle pulse
    elseif Tt1==0 && Tt2>0  
        t = linspace(0,Tt2,1000);
        for i=1:length(t)
            if t(i)/Tt2<1
                u(i) = 2*( (t(i)/Tt2)-0.5*(t(i)/Tt2).^2 ); %falling Triangle pulse
            else
                u(i)=0;
            end
        end
    elseif Tt1>0 && Tt2>0   %Tt2 is the time from the maximum of probability to 0
        t = linspace(0,Tt2+Tt1,1000);
        u=zeros(length(t),1);
        for i=1:length(t)
            if t(i)<= Tt1
                u(i) = (t(i)/Tt1).^2*(Tt1/(Tt1+Tt2));
            elseif  t(i)<= Tt1 + Tt2 
                u(i) = 2*(((t(i)-Tt1)/Tt2)-0.5*((t(i)-Tt1)/Tt2).^2)*(Tt2/(Tt1+Tt2))+Tt1/(Tt1+Tt2);
            else 
                u(i) = 0;
            end
        end

    elseif (Tt1>0 && Tt2>0) && Ts>0   %Tt2 is the time from the maximum of probability to 0
        t = linspace(0,Tt2+Tt1+Ts,1000);
        u=zeros(length(t),1);
        for i=1:length(t)
            if t(i)<= Tt1
                u(i) = (t(i)/Tt1).^2*(Tt1/(Tt1+Tt2+2*Ts));
            
            elseif  t(i)<= Tt1 + Ts
                u(i) = ((t(i)-Tt1)/Ts)*(2*Ts)/(2*Ts+Tt1+Tt2)+Tt1/(2*Ts+Tt1+Tt2);

            elseif  t(i)<= Tt1 + Tt2 +Ts
                u(i) = 2*(((t(i)-Tt1)/Tt2)-0.5*((t(i)-Tt1)/Tt2).^2)*(Tt2/(Tt1+Tt2+2*Ts))+(Tt1+2*Ts)/(Tt1+Tt2+2*Ts);
            else 
                u(i) = 0;
            end
        end

        elseif Tc>0
            t = linspace(0,Tc,1000);
            u=(-cos(pi*t/Tc)+1)/2;
    end
[u,index]=unique(u);
t=t(index);
crystal_par.t = t;
crystal_par.u = u;

% The parameters that will typically be the same for all SiPMs.





% The measured parameters of the SiPM at V_{ob} = V_{ob0}. Many a times,
% the values at V_{ob0} are the only ones that are available. Therefore,
% creating a separate structure for them.

% All the error bars correspond to 95 % confidence intervals, which corresponds to 2*sigma

%qe = [0.157 0.172]; % The internal microcell QE
qe = [pde 1];

%G_mean = [8.3*1e5 7.9*1e5]; % Average gain of the SiPM
%G_var = [0.5*1e5 0.5*1e5]/2; %Standard deviation of the gain

rdc_mean = [5.6*1e6 5.1*1e6]; % Rate of dark noise
rdc_var = [0.1e6 0.1e6]/2; %Standard deviation of dark noise
tau_phd1_mean = [15.3e-9 12.8e-9]; % The mean recovery time of the SiPM
tau_phd1_var = [1e-9 1e-9]/2; %Standard dev. of the recovery time

%Vob_mean = [ 1.26 1.39]; % The measured voltage over breakdown V_{ob0}
%Vob_var = [0.05 0.05]/2; %Symbols as defined in the paper

Vob_mean = [V 5.];
Vob_var = [0. 0.];

%Rq_mean = [145.3e3 140.4e3];
%Rq_var = [ 0.5e3 0.5e3]/2;


Rq_mean = [T_rec/(1.7e6*1.6e-19/3) 0];
Rq_var = [ 0 0]/2;

Rd_mean = [1e3 1e3];
Rd_var = [0 0];

taumr_mean = [1e-12 1e-12];
taumr_var = [0 0];

%%%
%Cd_mean = [77e-15 0];
%Cq_mean = [23e-15 0];

Cd_mean = [1.7e6*1.6e-19/3 0];
Cq_mean = [0 0];

Cg_mean = [28e-12 0];

Cd_var = [0 0];
Cq_var = [0 0];
Cg_var = [0 0];
%%%

V0_mean=[0,0];%Trigger variables
V0_var=[0,0];

%Vch_mean=[2.68,2];
Vch_mean=[2.68,2];

Vch_var=[0,0];

%Noise Parameters

%K_noise_mean=[10,0];
%K_noise_var=[0,0];
%Alpha_noise_mean=[1,0];
%Alpha_noise_var=[0,0];


if flags.oct==1
    P_ct=Prob;
else 
    P_ct=0;
end

if flags.ap==1
    P_ap=Prob;
else 
    P_ap=0;
end

if flags.doct==1
    P_dcr=Prob;
else 
    P_dcr=0;
end


P_prompt_ct_mean=[P_ct, 0.];
P_prompt_ct_var=[0., 0.];

P_ap_mean=[P_ap, 0.];
P_ap_var=[0., 0.];

P_delayed_ct_mean=[P_dcr, 0.];
P_delayed_ct_var=[0., 0.];

%Tbulk_mean=[10e-9,10e-9];
Tbulk_mean=[20e-9,10e-9];

Tbulk_var=[0,0];
oct_delay_mean=[1e-12,3e-9]; % 1 ps
oct_delay_var=[0,0];
ap_delay_mean=[3e-9,0];
ap_delay_var=[0,0];
doct_delay_mean=[3e-9,0];
doct_delay_var=[0,0];


Ratio_short_long_ap_mean=[0.1,0];
Ratio_short_long_ap_var=[0,0];

f = 1; % Collection efficiency of the SiPM
%f = 0.138;
%f = Prob;

meas_sipm_par(sipm_no) = struct('qe',0,'Nbar_oct_mean',0, 'Nbar_oct_var',0,...
    'Nbar_ap_short_mean',0,'Nbar_ap_short_var',0,'Nbar_ap_long_mean',0,'Nbar_ap_long_var',0,...
    'G_mean',0,'G_var',0,'rdc_mean',0,'rdc_var',0,...
    'tau_phd1_mean',0, 'tau_phd1_var',0,'tau_ap_short_mean',0, 'tau_ap_short_var',0,...
    'tau_ap_long_mean',0, 'tau_ap_long_var',0,...
    'V_ob_mean',0,'V_ob_var',0,'f',0,'Rq_mean',0,'Rq_var',0,...
    'pde_ratio',0,'ov_ratio_pde',0,'V0_mean',0,'V0_var',0,'Vch_mean',0,'Vch_var',0,'K_noise_mean',0,'K_noise_var',0,...
    'Alpha_noise_mean',0,'Alpha_noise_var',0,'Cd_mean',0,'Cq_mean',0,'Cg_mean',0,'Cd_var',0,'Cq_var',0,'Cg_var',0,...
    'Rd_mean',0,'Rd_var',0,'Tbulk_mean',0,'Tbulk_var',0,'oct_delay_mean',0,'oct_delay_var',0,'ap_delay_mean',0,...
    'ap_delay_var',0,'dct_delay_mean',0,'dct_delay_var',0,'Ratio_short_long_ap_mean',0,'Ratio_short_long_ap_var',0);
%     'xmax',1.5e-3,'xmin',-1.5e-3,'ymax',1.5e-3,'ymin',-1.5e-3);

% VM:

% qe:
% Nbar_oct_mean:
% Nbar_oct_var:
% Nbar_ap_short_mean:
% Nbar_ap_short_var:
% Nbar_ap_long_mean:
% Nbar_ap_long_var:
% G_mean:
% G_var:
% rdc_mean:
% rdc_var:
% tau_phd1_mean:
% tau_phd1_var:
% tau_ap_short_mean:
% tau_ap_short_var:
% tau_ap_long_mean:
% tau_ap_long_var:
% V_ob_mean:
% V_ob_var:
% f:
% Rq_mean:
% Rq_var:
% pde_ratio:
% ov_ratio_pde:
% V0_mean:
% V0_var:
% Vch_mean:
% Vch_var:
% K_noise_mean:
% K_noise_var:
% Alpha_noise_mean:
% Alpha_noise_var:
% Cd_mean:
% Cq_mean:
% Cg_mean:
% Cd_var:
% Cq_var:
% Cg_var:
% Rd_mean:
% Rd_var:
% Tbulk_mean:
% Tbulk_var:
% oct_delay_mean:
% oct_delay_var:
% ap_delay_mean:
% ap_delay_var:
% doct_delay_mean:
% doct_delay_var:


%Ncells = 2704; % The number of cells in the SiPM
Ncells = 100;
%Ncells = 676;
%Ncells=100;
%tmeas = 500e-9; %The entire measurement time

%We now define the geometry of the system. The SiPMs are placed in the x-y
%plane.
xmax = 1.5e-3; %The maximum coordinate along the x axis
xmin = -1.5e-3; %The mimimum coordinate along the x axis
ymax = 1.5e-3; %Same for y axis.
ymin = -1.5e-3;

sipm_pix_xno = sqrt(Ncells);
sipm_pix_yno = sqrt(Ncells);
wid_mcell_x = (xmax - xmin)/sipm_pix_xno;
wid_mcell_y = (ymax - ymin)/sipm_pix_yno;

tau_aval=0.1e-9;
tau_phd2=0.1e-12;
tau_phr=0.1e-9;
%original value
Rs=15;

%Rs=500;



fixed_sipm_par = struct('Ncells',Ncells,'tau_aval',tau_aval,'tau_phd2',tau_phd2,'tau_phr',tau_phr,'Rs',Rs,...
    'wid_mcell_x',wid_mcell_x,'wid_mcell_y',wid_mcell_y,'sipm_pix_xno',sipm_pix_xno,...
    'sipm_pix_yno',sipm_pix_xno,'tmeas',tmeas);

% VM:
% Ncells: Number of cells in the SiPM
% tau_aval: duration time of a charge avalanche 
% tau_phd2: Unused parameter from original code
% tau_phr: Unused parameter from original code
% Rs: load resistance of the readout circuit
% wid_mcell_x: width of single microcell in the x axis
% wid_mcell_y: width of single microcell in the y axis
% sipm_pix_xno: total number of microcells in the x axis
% sipm_pix_yno: total number of microcells in the y axis
% tmeas: Time of measurement


% Here we assume that all SiPM measurement parameters are the same, which is typically true.
% Thus, we just put in the measured values for the first SiPM and at the
% completion of the for loop, provide the same values for all other SiPMs.
% However, if that is not the case, simply change the for loop index from 1
% to number of SiPMs.

for sipm_ind = 1:1
    
    %     input_tmp = input('Please enter the input temperature at which the SiPM is operating: ','s');
    %     res_tmp = input('Do you have any measurements of SiPM parameter variation with temperature? Please type yes or no: ','s');
    %     if(strcmp(res_tmp,'yes'))
    %         %Please add the temperature varying component code here.
    %     end
    
    %     display('Reading the SiPM parameters at V_{ob} = V_{ob0} from file sipm_par_Vob0.dat')
    
    %Some initializations for variations with overvoltage
    
    pde_ratio = 1;
    ov_ratio_pde =1;
    ov_ratio_dc = 1;
    dc_ratio = 1;
    oct_ratio =1;
    ov_ratio_oct=1;
    ap_ratio=1;
    ov_ratio_ap=1;
    
    if(flag_default_variation_parameters_overvoltage)
        flags.overvoltage_pde = 1; 
        flags.overvoltage_pde_lin = 1; %PDE varies with overvoltage linearly
        flags.overvoltage_dc = 1; 
        flags.overvoltage_dc_lin = 1; % Dark current rate varies with overvoltage linearly
        flags.overvoltage_oct = 1;
        flags.overvoltage_oct_quad = 1; % Optical crosstalk varies with overvoltage quadratically
        flags.overvoltage_ap = 1;
        flags.overvoltage_ap_quad = 1; %Afterpulsing varies with overvoltage quadratically
    else
        
        
        display('---------------------------');
        display(['           SiPM ' num2str(sipm_ind)]);
        display('---------------------------');
        display('Detector PDE variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for PDE measurements with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_pde = 1;
            flags.overvoltage_pde_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: '); %Enter the value of PDE at V_{ob0} as the last entry
            %         no_meas = 5;
            filename = ['measurements\pde_ov_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen(filename,'r');
            if(fid == -1)
                error('File for PDE measurements vs overvoltage not found');
            end
            ov = fscanf(fid,'%f',no_meas); % The overvoltages written  as V_{ob}/V_{obz}.
            pde = fscanf(fid,'%f',no_meas); % The pde written as pde(V_{ob})/V_{obz}.
            ov_ratio_pde = ov/ov(no_meas);
            pde_ratio = pde/pde(no_meas);
            fclose(fid);
        else
            res = input('Would you like to use the usual linear relationship between overvoltage and PDE (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_pde = 1;
                flags.overvoltage_pde_lin = 1;
            else
                flags.overvoltage_pde = 0;
            end
        end
        
        display(' ');
        display('Dark current variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for dark current rate measurements with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_dc = 1;
            flags.overvoltage_dc_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: '); %Enter the value of PDE at V_{ob0} as the last entry
            filename = ['measurements\dc_ov_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen(filename,'r');
            if(fid == -1)
                error('File for dark current measurements vs overvoltage not found');
            end
            ov = fscanf(fid,'%f',no_meas); % The overvoltages written  as V_{ob}/V_{obz}.
            dc = fscanf(fid,'%f',no_meas); % The pde written as pde(V_{ob})/V_{obz}.
            ov_ratio_dc = ov/ov(no_meas);
            dc_ratio = dc/dc(no_meas);
            fclose(fid);
        else
            res = input('Would you like to use the usual linear relationship between overvoltage and dark current (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_dc = 1;
                flags.overvoltage_dc_lin = 1;
            else
                flags.overvoltage_dc = 0;
            end
        end
        
        display(' ');
        display('Optical Crosstalk probability variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for optical cross talk probability with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_oct = 1;
            flags.overvoltage_oct_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: ');
            filename = ['measurements\pde_oct_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen('filename','r');
            ov = fscanf(fid,'%f',no_meas);
            Nbar_oct = fscanf(fid,'%f',no_meas);
            fclose(fid);
            ov_ratio_oct = ov/ov(no_meas);
            oct_ratio = Nbar_oct/Nbar_oct(no_meas);
        else
            res = input('Would you like to use the usual quadratic relationship between overvoltage and after pulsing probability (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_oct = 1;
                flags.overvoltage_oct_quad = 1;
            else
                flags.overvoltage_oct = 0;
            end
        end
        
        display(' ');
        display('Afterpulsing probability variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for afterpulsing probability with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            
            flags.overvoltage_ap = 1;
            flags.overvoltage_ap_meas = 1;
            filename = ['measurements\pde_ap_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen('filename','r');
            ov = fscanf(fid,'%f',no_meas);
            Nbar_ap = fscanf(fid,'%f',no_meas);
            fclose(fid);
            ov_ratio_ap = ov/ov(no_meas);
            ap_ratio = Nbar_ap/Nbar_ap(no_meas);
        else
            res = input('Would you like to use the usual quadratic relationship between overvoltage and after pulsing probability (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_ap = 1;
                flags.overvoltage_ap_quad = 1;
            else
                flags.overvoltage_ap = 0;
            end
        end
        
        disp('');
        disp('');
    end
    
    meas_sipm_par(sipm_ind).qe = qe(sipm_ind);
    
    %meas_sipm_par(sipm_ind).G_mean = G_mean(sipm_ind);
    %meas_sipm_par(sipm_ind).G_var = G_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).rdc_mean = rdc_mean(sipm_ind);
    meas_sipm_par(sipm_ind).rdc_var = rdc_var(sipm_ind);
    meas_sipm_par(sipm_ind).tau_phd1_mean = tau_phd1_mean(sipm_ind);
    meas_sipm_par(sipm_ind).tau_phd1_var = tau_phd1_var(sipm_ind);
   
    meas_sipm_par(sipm_ind).V_ob_mean = Vob_mean(sipm_ind);
    meas_sipm_par(sipm_ind).V_ob_var = Vob_var(sipm_ind);
    meas_sipm_par(sipm_ind).Rq_mean = Rq_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Rq_var = Rq_var(sipm_ind);
    meas_sipm_par(sipm_ind).pde_ratio = pde_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_pde = ov_ratio_pde;
    meas_sipm_par(sipm_ind).dc_ratio = dc_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_dc = ov_ratio_dc;
    meas_sipm_par(sipm_ind).oct_ratio = oct_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_oct = ov_ratio_oct;
    meas_sipm_par(sipm_ind).ap_ratio = ap_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_ap = ov_ratio_ap;
    meas_sipm_par(sipm_ind).f = f;
    
    
    meas_sipm_par(sipm_ind).Cd_mean=Cd_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Cq_mean=Cq_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Cg_mean=Cg_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Rd_mean=Rd_mean(sipm_ind);
    meas_sipm_par(sipm_ind).taumr_mean=taumr_mean(sipm_ind);
    
    meas_sipm_par(sipm_ind).Cd_var=Cd_var(sipm_ind);
    meas_sipm_par(sipm_ind).Cq_var=Cq_var(sipm_ind);
    meas_sipm_par(sipm_ind).Cg_var=Cg_var(sipm_ind);
    meas_sipm_par(sipm_ind).Rd_var=Rd_var(sipm_ind);
    meas_sipm_par(sipm_ind).taumr_var=taumr_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).V0_mean = V0_mean(sipm_ind);%Trigger parameters
    meas_sipm_par(sipm_ind).V0_var = V0_var(sipm_ind);
    meas_sipm_par(sipm_ind).Vch_mean = Vch_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Vch_var = Vch_var(sipm_ind);
    
    %meas_sipm_par(sipm_ind).K_noise_mean=K_noise_mean(sipm_ind); %Noise parameters
    %meas_sipm_par(sipm_ind).K_noise_var=K_noise_var(sipm_ind);
    %meas_sipm_par(sipm_ind).Alpha_noise_mean=Alpha_noise_mean(sipm_ind);
    %meas_sipm_par(sipm_ind).Alpha_noise_var=Alpha_noise_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).P_prompt_ct_mean=P_prompt_ct_mean(sipm_ind);
    meas_sipm_par(sipm_ind).P_prompt_ct_var=P_prompt_ct_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).P_ap_mean=P_ap_mean(sipm_ind);
    meas_sipm_par(sipm_ind).P_ap_var=P_ap_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).P_delayed_ct_mean=P_delayed_ct_mean(sipm_ind);
    meas_sipm_par(sipm_ind).P_delayed_ct_var=P_delayed_ct_var(sipm_ind);
    
    meas_sipm_par(sipm_ind).Tbulk_mean=Tbulk_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Tbulk_var=Tbulk_var(sipm_ind);
    meas_sipm_par(sipm_ind).oct_delay_mean=oct_delay_mean(sipm_ind);
    meas_sipm_par(sipm_ind).oct_delay_var=oct_delay_var(sipm_ind);
    meas_sipm_par(sipm_ind).ap_delay_mean=ap_delay_mean(sipm_ind);
    meas_sipm_par(sipm_ind).ap_delay_var=ap_delay_var(sipm_ind);
    meas_sipm_par(sipm_ind).doct_delay_mean=doct_delay_mean(sipm_ind);
    meas_sipm_par(sipm_ind).doct_delay_var=doct_delay_var(sipm_ind);

    meas_sipm_par(sipm_ind).Ratio_short_long_ap_mean=Ratio_short_long_ap_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Ratio_short_long_ap_var=Ratio_short_long_ap_var(sipm_ind);
    
    
    
end

% Now let all the SiPMs have the same property that the first SiPM has
for sipm_ind = 2:sipm_no
    meas_sipm_par(sipm_ind) = meas_sipm_par(1);
end


disp('Check the values:')
if Td>0 && Tr==0 && Ts==0
    tipo='exponential';
    T=['Td: ' num2str(Td*10^9) ' ns'];
    T_cocient=['Td/Trec: ' num2str(Td/T_rec)];
end
if Td>0 && Tr>0 && Ts==0
    tipo='biexponential';
    T=['Td: ' num2str(Td*10^9) ' ns, Tr:' num2str(Tr*10^9) ' ns'];
    T_cocient=['Td/Trec: ' num2str(Td/T_rec) ' , Tr/Td: ' num2str(Td/Tr)];
end
if Td==0 && Tr==0 && Ts>0
    tipo='squared';
    T=['Ts: ' num2str(Ts*10^9) ' ns'];
    T_cocient=['Ts/Trec: ' num2str(Ts/T_rec)];
end
if Td==0 && Tr==0 && Ts==0 && Tt1>0 && Tt2==0
    tipo='rising triangular';
    T=['Tt1: ' num2str(Tt1*10^9) ' ns'];
    T_cocient=['Tt1/Trec: ' num2str(Tt1/T_rec)];
end
if Td==0 && Tr==0 && Ts==0 && Tt1==0 && Tt2>0
    tipo='decreasing triangular';
    T=['Tt2: ' num2str(Tt2*10^9) ' ns'];
    T_cocient=['Tt2/Trec: ' num2str(Tt2/T_rec)];
end
if Td==0 && Tr==0 && Ts==0 && Tt1>0 && Tt2>0
    tipo='triangular';
    T=['Tt1: ' num2str(Tt1*10^9) ' ns, Tt2: ' num2str(Tt2*10^9) ' ns'];
    T_cocient=['Tt1/Trec: ' num2str(Tt1/T_rec) ', Tt2/Trec: ' num2str(Tt2/T_rec) ' , Tt1/Tt2: ' num2str(Tt1/Tt2)];
end

if Ts>0 && Tt1>0 && Tt2>0
    tipo='Semi-rectangular';
    T=['Tt1: ' num2str(Tt1*10^9) ' ns, Ts: ' num2str(Ts*10^9) ' ns, Tt2: ' num2str(Tt2*10^9) ' ns'];
    T_cocient=['Tt1/Trec: ' num2str(Tt1/T_rec) ', Ts/Trec: ' num2str(Ts/T_rec) ', Ts/(Tt1+Tt2): ' num2str(Ts/(Tt1+Tt2)) ' , Tt1/Tt2: ' num2str(Tt1/Tt2)];
end

if Tc>0 
    tipo='Semi-circular';
    T=['Tc: ' num2str(Tc*10^9) ' ns'];
    T_cocient=['Tc/Trec: ' num2str(Tc/T_rec) ];
end
if Poisson ==1
    Poisson='True';
else
    Poisson='False';
end
if flags.simplified_readout ==1
    Simplified='True';
else
    Simplified='False';
end
disp(['Light pulse shape:' tipo])
disp([T ' '])
disp(['Poissonnian fluctuations: ' Poisson])
disp('SIPM PARAMETERSS')
disp(['Simplified read-out: ' Simplified ])
display(['Pixel number: ' num2str(Ncells)])
disp(['Yield: ' num2str(Qsc) ' photon/keV'])
disp(['Geometric efficiency: ' num2str(f)])
disp(['Recovery Time: ' num2str(T_rec*10^9) ' ns'])
disp(['PDE: ' num2str(pde)])
disp(['U: ' num2str(V) ' V'])
disp(['U0: ' num2str(V0_mean(1)) ' V'])
disp(['Uch: ' num2str(Vch_mean(1)) ' V'])
disp('NOISE PARAMETERS')
disp(['Pct: ' num2str(P_prompt_ct_mean(1)*100) ' %'])
disp(['Pap: ' num2str(P_ap_mean(1)*100) ' %'])
disp(['Pdct: ' num2str(P_delayed_ct_mean(1)*100) ' %'])
disp(['Tbulk: ' num2str(Tbulk_mean(1)*10^9) ' ns'])
disp('MODEL PARAMETERS')
disp([T_cocient ' '])
disp(['U/Uch: ' num2str(V/Vch_mean(1))])

wait=0;
if wait==1
    pause(45)
end

info = {
    'Light pulse shape:', tipo;
    'T:', T;
    'Poissonnian fluctuations:', Poisson;
    'SIPM PARAMETERS', '';
    'Simplified read-out:', Simplified;
    'Pixel number:', Ncells;
    'Yield (photon/keV):', Qsc;
    'Geometric efficiency:', f;
    'Recovery Time (ns):', T_rec*1e9;
    'PDE:', pde;
    'U (V):', V;
    'U0 (V):', V0_mean(1);
    'Uch (V):', Vch_mean(1);
    'NOISE PARAMETERS', '';
    'Pct (%):', P_prompt_ct_mean(1)*100;
    'Pap (%):', P_ap_mean(1)*100;
    'Pdct (%):', P_delayed_ct_mean(1)*100;
    'Tbulk (ns):', Tbulk_mean(1)*1e9;
    'MODEL PARAMETERS', '';
    'T_cocient:', T_cocient;
    'U/Uch:', V/Vch_mean(1);
    'Fecha:', datestr(now);
};
