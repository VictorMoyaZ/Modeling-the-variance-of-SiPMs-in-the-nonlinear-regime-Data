function [ph_fired,phot_pot] = sipm_sim_array(no_timesteps,Nph_t,meas_sipm_par,fixed_sipm_par, crystal_par, flags,s_mcell,s_mcell_hit,s_mcell_hit_first_gains  )
% The code to simulate the SiPM operation.
% The following are the input arguments:
% no_timesteps: The number of time steps into which you would like to
% divide the total measurement time
% Nph_t: The timestamps of the scintillation photons
% meas_sipm_par: The measured parameters of the SiPM. Please see
% read_input_array.m
% fixed_sipm_par: The fixed SiPM parameters. Again see read_input_array.m
% crystal_par: The crsytal decay and rise times. Please see read_input_array.m
% flags: The flags to turn on and off the different SiPM mechanisms.
% Default value is 1. Please see read_input_array.m for more details.


% Read in all the input variables to the local variables in this code

%%% VM: First the measured parameters of the SiPM and its stochastic generation
%
qe = meas_sipm_par.qe;
charge_elec = 1.6e-19;

%G_mean = meas_sipm_par.G_mean;
%G_var = meas_sipm_par.G_var;


rdc_mean = meas_sipm_par.rdc_mean;
rdc_var = meas_sipm_par.rdc_var;
V_ob_mean = meas_sipm_par.V_ob_mean;
V_ob_var = meas_sipm_par.V_ob_var;
Rq_mean = meas_sipm_par.Rq_mean;
Rq_var = meas_sipm_par.Rq_var;
f = meas_sipm_par.f;

V0_mean = meas_sipm_par.V0_mean ;%VM: Trigger parameters
V0_var = meas_sipm_par.V0_var ;
Vch_mean = meas_sipm_par.Vch_mean;
Vch_var = meas_sipm_par.Vch_var;




%G = G_mean + randn(1)*G_var;
rdc = rdc_mean + randn(1)*rdc_var;

V_ob = V_ob_mean + randn(1)*V_ob_var;
Rq = Rq_mean + randn(1)*Rq_var;

%VM: Trigger parameters

V0=V0_mean + randn(1)*V0_var;
Vch=Vch_mean + randn(1)*Vch_var;

Norm_trigger=1-exp(-(V_ob-V0)/Vch);
%%%


% The overvoltage variations have been deleted

    
%%% Now the fixed parameters of the SiPM
Ncells = fixed_sipm_par.Ncells;
Rs = fixed_sipm_par.Rs;
sipm_pix_xno = fixed_sipm_par.sipm_pix_xno;
sipm_pix_yno = fixed_sipm_par.sipm_pix_yno;
tmeas = fixed_sipm_par.tmeas;
%%%

%%% Crystal parameters
Td = crystal_par.Td;
Qsc = crystal_par.Qsc;

%%%


%%%VM: Electric response parameters

if flags.simplified_readout
    
    Cd =  meas_sipm_par.Cd_mean + randn(1)* meas_sipm_par.Cd_var;
    Cq = meas_sipm_par.Cq_mean + randn(1)* meas_sipm_par.Cq_var;

    G=(Cd+Cq)*V_ob/charge_elec;


    tau1 = charge_elec*G*Rq/V_ob; % VM: revisar
    tau2 = 1;
    tau3 = 1;
    tau4 = 1;

    taumr = 1;
    taumd = 1;
    tau_aval = 1;

    taucd1 = 1;
    taucd2 = 1;

    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0;

    b1 = 1/tau1;
    b2 = 0;
    b3 = 0;
    b4 = 0;
    b5 = 0;
    b6 = 0;
    b7 = 0;
    b8 = 0;
else
    Cd =  meas_sipm_par.Cd_mean + randn(1)* meas_sipm_par.Cd_var;
    Cq = meas_sipm_par.Cq_mean + randn(1)* meas_sipm_par.Cq_var;
    Cg = meas_sipm_par.Cg_mean + randn(1)* meas_sipm_par.Cg_var;
    
    G=(Cd+Cq)*V_ob/charge_elec;
    
    Rq = Rq_mean+ randn(1)*Rq_var;
    Rd = meas_sipm_par.Rd_mean + randn(1)* meas_sipm_par.Rd_var;
    Rqp = Rq/(Ncells-1);
    Cqp = Cq*(Ncells-1);
    Cdp = Cd*(Ncells-1);


    tau1 = charge_elec * G*Rq/V_ob;
    tau2 = Rq*Cq;
    tau3 = Rs*Ncells*Cd;
    tau4 = Rs*Cg;

    taumr = meas_sipm_par.taumr_mean + randn(1)* meas_sipm_par.taumr_var;
    taumd = Cd*(Rd*Rq/(Rd+Rq));
    tau_aval = taumd;

    taucd1 = ((tau1 + tau3 + tau4 ) + sqrt( (tau1 + tau3 + tau4 )^2 - 4 * ( tau1 * tau4 + tau2 * tau3))) / 2;
    taucd2 = ((tau1 + tau3 + tau4 ) - sqrt( (tau1 + tau3 + tau4 )^2 - 4 * ( tau1 * tau4 + tau2 * tau3))) / 2;

    a1 = taucd1*(taucd1 - tau2)/((taucd1 - taucd2) * ( taucd1 - taumd) * ( taucd1 - taumr));
    a2 = taucd2*(taucd2 - tau2)/((taucd2 - taucd1) * ( taucd2 - taumd) * ( taucd2 - taumr));
    a3 = taumd*(taumd - tau2)/((taumd - taucd1) * ( taumd - taucd2) * ( taumd - taumr));
    a4 = taumr*(taumr - tau2)/((taumr - taucd1) * ( taumr - taucd2) * ( taumr - taumd));

    b1 = tau1/( (tau1 - taumr)*(tau1 - taumd));
    b2 = taumd/( (taumd - tau1)*(taumd - taumr));
    b3 = taumr/( (taumr - tau1)*(taumr - taumd));
    b4 = (taucd1 - tau2)*taucd1/( (taucd1 - taucd2)*(taucd1 - tau1)*(taucd1 - taumd)*(taucd1 - taumr));
    b5 = (taucd2 - tau2)*taucd2/( (taucd2 - taucd1)*(taucd2 - tau1)*(taucd2 - taumd)*(taucd2 - taumr));
    b6 = (taumd - tau2)*taumd/( (taumd - taucd1)*(taumd - tau1)*(taumd - taucd2)*(taumd - taumr));
    b7 = (taumr - tau2)*taumr/( (taumr - taucd1)*(taumr - tau1)*(taumr - taumd)*(taumr - taucd2));
    b8 = (tau1 - tau2)*tau1/( (tau1 - taucd1)*(tau1 - taucd2)*(tau1 - taumd)*(tau1 - taumr)); 
end
%%%
%%% VM: Noise parameters 


oct_delay = meas_sipm_par.oct_delay_mean + randn(1)*meas_sipm_par.oct_delay_var;
ap_delay = meas_sipm_par.ap_delay_mean + randn(1)*meas_sipm_par.ap_delay_var;
doct_delay = meas_sipm_par.doct_delay_mean + randn(1)*meas_sipm_par.doct_delay_var; 

%K_noise_mean=meas_sipm_par.K_noise_mean; %Noise parameters
%K_noise_var=meas_sipm_par.K_noise_var;
%Alpha_noise_mean=meas_sipm_par.Alpha_noise_mean;
%Alpha_noise_var=meas_sipm_par.Alpha_noise_var;

%K_noise=K_noise_mean + randn(1)*K_noise_var;
%Alpha_noise=Alpha_noise_mean+randn(1)*Alpha_noise_var;

P_prompt_ct_mean=meas_sipm_par.P_prompt_ct_mean;
P_prompt_ct_var=meas_sipm_par.P_prompt_ct_var;
P_prompt_ct=P_prompt_ct_mean+randn(1)*P_prompt_ct_var;
lambda_prompt_ct_operative=-log(1-P_prompt_ct);

P_ap_mean=meas_sipm_par.P_ap_mean;
P_ap_var=meas_sipm_par.P_ap_var;
P_ap_operative=P_ap_mean+randn(1)*P_ap_var;    

P_delayed_ct_mean=meas_sipm_par.P_delayed_ct_mean;
P_delayed_ct_var=meas_sipm_par.P_delayed_ct_var;
P_delayed_ct=P_delayed_ct_mean+randn(1)*P_delayed_ct_var;
lambda_delayed_ct_operative=-log(1-P_delayed_ct);

if flags.Biexponential_ap_distibution
    lambda_ap_operative=-log(1-P_ap_operative);

    Ratio_short_long_ap=meas_sipm_par.Ratio_short_long_ap_mean+randn(1)*meas_sipm_par.Ratio_short_long_ap_var;

else

    Tbulk=meas_sipm_par.Tbulk_mean+randn(1)*meas_sipm_par.Tbulk_var;

    [P_ap,T_ap,I_ap] = distribution_integrator(-1,Tbulk,30,ap_delay);
    [P_dct,T_dct,I_dct] = distribution_integrator(-0.5,Tbulk,30,doct_delay);

    Afterpulsing_Integrated_distribution_1=I_ap;
    Afterpulsing_Integrated_distribution_2=Afterpulsing_distribution_integrator(-1,Tbulk,tau1,ap_delay,V_ob,V0,Vch);

    lambda_ap_operative=-Afterpulsing_Integrated_distribution_1*log(1-P_ap_operative)/Afterpulsing_Integrated_distribution_2;
end

%%%



% VM:recovery variables have been replaced by trigger








%% Running the simulation
%%% VM: Tsample calculation

if (flags.output_voltage == 1 && no_timesteps > 1)
    num_time_samples = no_timesteps;
    tmp = linspace(0,1 - exp(-6),num_time_samples);
    Tsamp = -Td * log(1 - tmp);
    if(Tsamp(num_time_samples) > tmeas)
        Tsamp(num_time_samples)= tmeas;
    end
    %    Tsamp = linspace(0,tmeas,num_time_samples);
    delt = Tsamp(2:num_time_samples) - Tsamp(1:num_time_samples-1);
    delt(num_time_samples) = tmeas - Tsamp(num_time_samples);
else
    num_time_samples = 1; %ceil(tmeas/tau_phd1);
    Tsamp = 0;
    delt = tmeas;
end

%%%

%%% VM: Creation of data structures
%mcell(sipm_pix_xno,sipm_pix_yno,num_time_samples) = struct('count',0,'t',[],'gain',[],'flag',[]);
%mcell_hit(sipm_pix_xno,sipm_pix_yno) = struct('count',0,'t',[],'flag',[]);
%mcell_hit_first_gains(sipm_pix_xno,sipm_pix_yno) = struct('first_photon',0,'noise_count',1,'first_gain',0,'first_noise_time',[],'first_noise_origin_i',[],'first_noise_origin_j',[]);
%phot_pot = struct('count',0,'t',[], 'flag',[],'i',[],'j',[],'gain',[],'analysed_count',0,'first_gain_mean',0);


mcell = repmat(s_mcell, sipm_pix_xno, sipm_pix_yno, num_time_samples);
mcell_hit = repmat(s_mcell_hit, sipm_pix_xno, sipm_pix_yno);
mcell_hit_first_gains = repmat(s_mcell_hit_first_gains, sipm_pix_xno, sipm_pix_yno);
phot_pot = struct('count',0,'t',[], 'flag',[],'i',[],'j',[],'gain',[],'mu',[],'sigma',[],'analysed_count',0,'first_gain_mean',0,'trigger',[]);


for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        mcell_hit_first_gains(i,j).first_photon=0;
        mcell_hit_first_gains(i,j).noise_count=1;
        mcell_hit_first_gains(i,j).first_gain=0;
        mcell_hit_first_gains(i,j).first_noise_time=[];
        mcell_hit_first_gains(i,j).first_noise_origin_i=[];
        mcell_hit_first_gains(i,j).first_noise_origin_j=[];
    end
end


%%%

for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        mcell_hit(i,j).count = 0;
        mcell_hit(i,j).t = [];
        for k = 1:num_time_samples,
            mcell(i,j,k).count = 0;
            mcell(i,j,k).t = [];
        end
    end
end



% Now assume that when the photons travel to the SiPM, the order between them is maintained.



% Of the generated photons, some photons reach the SiPM and some do not.

Nph = length(Nph_t);
coll_success = rand(Nph,1);
ind = (coll_success < f);
Nph_sipm = Nph_t(ind);
Nph_sipm_no = length(Nph_sipm);
Nph_sipm = sort(Nph_sipm,'ascend');

%
phot_rej = 0;
lost_phot = 0;

%%% VM:Assigning each photon its arrival time and cell

for ph_ind = 1:Nph_sipm_no,
    i = ceil(rand(1)*sipm_pix_xno);
    j = ceil(rand(1)*sipm_pix_yno);
    
    if(rand(1) < qe)
        if (mcell_hit(i,j).count > 0 && flags.deadtime )
            if (Nph_sipm(ph_ind) - mcell_hit(i,j).t(mcell_hit(i,j).count) > tau_aval)
                %if mcell_hit(i, j).count>0
                %    break
                %end
                mcell_hit(i, j).t = [mcell_hit(i, j).t Nph_sipm(ph_ind)];
                mcell_hit(i, j).count = mcell_hit(i, j).count + 1;
                mcell_hit(i,j).flag = [mcell_hit(i, j).flag 1];
                phot_pot.t = [phot_pot.t Nph_sipm(ph_ind)];
                phot_pot.flag = [phot_pot.flag 1];
                phot_pot.i = [phot_pot.i i];
                phot_pot.j = [phot_pot.j j];
                phot_pot.gain = [phot_pot.gain 0];
                phot_pot.count = phot_pot.count + 1;

                %% VM
                if isempty(mcell_hit(i,j).first_time)
                    mcell_hit(i,j).first_time=[mcell_hit(i,j).first_time Nph_sipm(ph_ind)];
                end

            else
                phot_rej = phot_rej + 1;
            end
        elseif (~flags.deadtime || mcell_hit(i,j).count == 0)

%             if mcell_hit(i, j).count>0
%                continue
%             end
            mcell_hit(i, j).t = [mcell_hit(i, j).t Nph_sipm(ph_ind)];
            mcell_hit(i, j).count = mcell_hit(i, j).count + 1;
            phot_pot.t = [phot_pot.t Nph_sipm(ph_ind)];
            phot_pot.flag = [phot_pot.flag 1];
            phot_pot.i = [phot_pot.i i];
            phot_pot.j = [phot_pot.j j];
            phot_pot.gain = [phot_pot.gain 0];
            phot_pot.mu = [phot_pot.mu 0];
            phot_pot.sigma = [phot_pot.sigma 0];
            phot_pot.count = phot_pot.count + 1;
            
        end
    end
end


%%%


%%% VM:Generating Darck Current events

if(flags.dc)
    % The times at which the dark current will trigger avalanches
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            prob_darkcurr = rdc*tmeas/Ncells;
            Ndc = poissrnd(prob_darkcurr); %Stochastic
            %Ndc = uint8(prob_darkcurr);
            if(Ndc > 0)
                darkcurr_time = rand(Ndc,1)*tmeas;
                mcell_hit(i,j).t = [mcell_hit(i, j).t darkcurr_time'];
                mcell_hit(i,j).count = mcell_hit(i, j).count + Ndc;
                phot_pot.t = [phot_pot.t darkcurr_time'];
                phot_pot.flag = [phot_pot.flag 2*ones(Ndc,1)'];
                phot_pot.i = [phot_pot.i i*ones(Ndc,1)'];
                phot_pot.j = [phot_pot.j j*ones(Ndc,1)'];
                phot_pot.gain = [phot_pot.gain 0*ones(Ndc,1)'];
                phot_pot.count = phot_pot.count + Ndc;
            end
        end
    end
end


%%%

%%% VM:temporarily sort of events
for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        [mcell_hit(i,j).t] = sort(mcell_hit(i,j).t,'ascend');
        %             mcell_hit(i,j).flag = mcell_hit(i,j).flag(order);
    end
end


[phot_pot.t,order] = sort(phot_pot.t, 'ascend');
phot_pot.flag = phot_pot.flag(order);
phot_pot.i = phot_pot.i(order);
phot_pot.j = phot_pot.j(order);
phot_pot.gain = phot_pot.gain(order);


%%%

% Now that we know which photons are hitting which microcell, and since
% the various parameters of the SiPM are themselves function of
% over-voltage which in turn is a function of output voltage, therefore, we need to sample our time interval so that
% at the end of different time intervals, we evaluate the  output voltage
% and use that to redetermine the various SiPM parameters.


%%% VM: Simulation loop of events

Vs = 0;
ind_t = 1;
phot_pot_samp_count = zeros(num_time_samples,1);
phot_glob_count = 0;
%     ap_count_2 = 0;
%     oct_count_2 = 0;
I_Rs = 0;

for t_samp = Tsamp, % VM: this first loop works on the time samples of simulation
    
    %%% VM: Calculus of Vs voltage at event time
    if(flags.output_voltage && ind_t > 1 && flags.simplified_readout == 0)
        Vs = (Rs) * I_Rs ; %(charge_time_int(ind_t-1)) * G * charge_elec/delt(ind_t-1);
        Vs_ind(ind_t) = Vs;
    end
    %%%
    
    %%% VM: Adding a new event to mcell and sorting it by its generation time
    
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            for k = 1:mcell_hit(i,j).count,
                if(mcell_hit(i,j).t(k) >= t_samp && mcell_hit(i,j).t(k) < t_samp + delt(ind_t)),
                    mcell(i,j,ind_t).t = [mcell(i, j,ind_t).t mcell_hit(i,j).t(k)];
                    mcell(i,j,ind_t).count = mcell(i, j,ind_t).count + 1;
                    mcell(i,j,ind_t).gain = [mcell(i, j,ind_t).gain 0];
                    %                         mcell(i,j,ind_t).flag = [mcell(i, j,ind_t).flag mcell_hit(i,j).flag];
                    phot_pot_samp_count(ind_t) = phot_pot_samp_count(ind_t) + 1;
                end
            end
        end
    end
    
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            
            [mcell(i,j,ind_t).t, order] = sort(mcell(i,j,ind_t).t,'ascend');
        end
    end
    
    
    %%%
    
    
    %%% VM:this second loop works on each event on the current time sample
    gain = 1;
    phot_count = 0;
    
    if(phot_pot_samp_count(ind_t) > 0) 
        photon_present = 1;
    else
        photon_present = 0;
    end
    
    while photon_present,  
        
        %VM: updating the event counter
        phot_count = phot_count + 1;
        phot_glob_count = phot_glob_count + 1;
        %
        %VM: obtaining cell's index
        i = phot_pot.i(phot_glob_count); 
        j = phot_pot.j(phot_glob_count);
        %
        
        %VM: getting the event time index in the current time stamp 
        for k = 1:mcell(i,j,ind_t).count,
            if( mcell(i,j,ind_t).t(k) == phot_pot.t(phot_glob_count))
                break;
            end
        end
        
        if(isempty(k))
            display('break here due to empty k');
        end
        %
        
        % VM: gain calculation
        if(flags.partial_recovery || flags.over_voltage_effect)
            
            
            I_Rq = 0; % Current through quenching resistance due to all the previous avalanches
            if(ind_t >1)
                for prev_ind_t = 1:ind_t-1
                    for prev_phot_ind = 1:mcell(i,j,prev_ind_t).count
                        prev_trigger_time = mcell(i,j,prev_ind_t).t(prev_phot_ind);
                        tp = mcell(i,j,ind_t).t(k) - prev_trigger_time;
                        I_Rq = I_Rq + (G*charge_elec) * mcell(i,j,prev_ind_t).gain(prev_phot_ind)* ...
                            ((b1 * exp(-tp/tau1) + b2*exp(-tp/taumd) + b3*exp(-tp/taumr)) + ...
                            (tau3/Ncells)*(b4*exp(-tp/taucd1) + b5*exp(-tp/taucd2) + b6*exp(-tp/taumd) + b7*exp(-tp/taumr)+ b8*exp(-tp/tau1)));
                    end
                end
            end
            if( k >1)
                for prev_phot_ind = 1:k-1
                    prev_trigger_time = mcell(i,j,ind_t).t(prev_phot_ind);
                    tp = mcell(i,j,ind_t).t(k) - prev_trigger_time;
                    I_Rq = I_Rq + (G*charge_elec) * mcell(i,j,ind_t).gain(prev_phot_ind)* ...
                        ((b1 * exp(-tp/tau1) + b2*exp(-tp/taumd) + b3*exp(-tp/taumr)) - ...
                        (tau3/Ncells)*(b4*exp(-tp/taucd1) + b5*exp(-tp/taucd2) + b6*exp(-tp/taumd) + b7*exp(-tp/taumr) + b8*exp(-tp/tau1)));
                    
                end
            end
            
            
            gain = (1 - Vs/V_ob - I_Rq*Rq/V_ob);
            if(gain < 0)
                gain = 0;
                lost_phot = lost_phot+1 ;
            end
            if(isnan(gain))
                display('gain is nan');
            end
            if(gain > 1)
                gain = 0;
            end
        else
            gain = 1;
        end
        
        if(isnan(gain))
            display('gain is nan');
        end
        %
        
        % Poisson fluctuation
       
        mu = gain ;       % media
        sigma = (0.6546*exp(-0.4424*gain*5.2)+0.0796)*gain; 
        
        AAA=randn(1);% desviación estándar
        gain = mu + sigma * AAA;
        
        

        % VM: trigger calculation 
        
        if (flags.trigger)
            %Trigger=(1-exp(-(gain*V_ob-V0)/Vch))/Norm_trigger;
            Trigger=(1-exp(-(mu*V_ob-V0)/Vch))/Norm_trigger;
            
        else
            Trigger=1;
        end
        
        %
        
        %VM: overvoltage recovery flags has been deleted
        
        %%% This adds noise events to the event data list
        
        % Noise Simulation
        BBB=rand(1);
        if( BBB < Trigger)
            % Check if the photon or dark current or
            % afterpulses ( due to previous photon triggered avalanches) or OCT triggers an avalanche. If so, only
            % then is there scope for OCT, afterpulse etc.
            % Afterpulses do not require the gain parameter
            % to be evaluated for a different overvoltage. We have
            % already considered that when evaluating prob(N_ap).
            
            phot_pot.gain(phot_glob_count) = gain;
            % phot_pot.mu(phot_glob_count) = BBB;
            % phot_pot.sigma(phot_glob_count) = AAA;
            % phot_pot.trigger(phot_glob_count) = Trigger;
            mcell(i,j,ind_t).gain(k) = mu;
            %mcell(i,j,ind_t).gain(k) = gain;
            
            
            
            
            % The times at which the after pulses occur. The afterpulses are
            % inhomogeneous Poisson point processes with intensity function being
            % an exponential with decay time tap, and the mean of the Poisson being
            % N_ap. Therefore, this event will occur only
            % after the avalanche that triggers it, and will be
            % following it in the list.
            
            % Afterpulsing simulation
            if(flags.ap)
                
                %prob=1-exp(-(V_ob*gain/K_noise)^(1+Alpha_noise));
                %lambda=-log(1-prob);
                %C=lambda/Afterpulsing_Integrated_distribution_1;
                %lamdba_2=C*Afterpulsing_Integrated_distribution_2;
                Nap = poissrnd(lambda_ap_operative*gain);
                
                
                
                if(Nap > 0)
                    ap_time = get_delayed_noise_time(Nap,P_ap,T_ap); %exprnd(tau_ap,Nap,1);
                    for kp = 1:Nap,
                        ap_time_ind = 0;
                        time_ap = ap_time(kp) + mcell(i,j,ind_t).t(k);
                        for tmp_ind_t = ind_t+1:num_time_samples,
                            if(time_ap < Tsamp(tmp_ind_t))
                                ap_time_ind = tmp_ind_t-1;
                                break;
                            end
                        end
                        if(time_ap > Tsamp(num_time_samples) && time_ap < tmeas)
                            ap_time_ind = num_time_samples;
                        end
                        if(ap_time_ind)
                            mcell(i,j,ap_time_ind).t = [mcell(i,j,ap_time_ind).t ap_time(kp) + mcell(i,j,ind_t).t(k)];
                            mcell(i,j,ap_time_ind).count = mcell(i,j,ap_time_ind).count + 1;
                            mcell(i,j,ap_time_ind).gain = [mcell(i,j,ap_time_ind).gain 0];
                            %                                 mcell(i,j,ap_time_ind).flag = [mcell(i,j,ap_time_ind).flag 3];
                            %ap_count_2 = ap_count_2 + 1;
                            % Re-sort to get the afterpulses at the correct place in
                            % the list.
                            [mcell(i,j,ap_time_ind).t, order] = sort(mcell(i,j,ap_time_ind).t,'ascend');
                            phot_pot.t = [phot_pot.t time_ap];
                            phot_pot.flag = [phot_pot.flag 3];
                            phot_pot.i = [phot_pot.i i];
                            phot_pot.j = [phot_pot.j j];
                            phot_pot.gain = [phot_pot.gain 0];
                            phot_pot.count = phot_pot.count + 1;
                            
                            phot_pot_samp_count(ap_time_ind) = phot_pot_samp_count(ap_time_ind) + 1;
                            
                            if mcell_hit_first_gains(i,j).first_photon==0 && (phot_pot.flag(phot_glob_count)==1 || mcell_hit_first_gains(i,j).noise_count<=length(mcell_hit_first_gains(i,j).first_noise_time))
                            ip=i;
                            jp=j;
                            mcell_hit_first_gains(ip,jp).first_noise_origin_i=[mcell_hit_first_gains(ip,jp).first_noise_origin_i i];
                            mcell_hit_first_gains(ip,jp).first_noise_origin_j=[mcell_hit_first_gains(ip,jp).first_noise_origin_j j];
                            mcell_hit_first_gains(ip,jp).first_noise_time=[mcell_hit_first_gains(ip,jp).first_noise_time time_ap];
                            

                            [mcell_hit_first_gains(ip,jp).first_noise_time,order] = sort(mcell_hit_first_gains(ip,jp).first_noise_time, 'ascend');
                            mcell_hit_first_gains(ip,jp).first_noise_origin_i = mcell_hit_first_gains(ip,jp).first_noise_origin_i(order);
                            mcell_hit_first_gains(ip,jp).first_noise_origin_j = mcell_hit_first_gains(ip,jp).first_noise_origin_j(order);
       

                            
                            end
                        end
                    end
                end
            end
            %
            
            % The OCT events occur within 1 ns in the neighboring 4
            % cells ( with uniform probability).
            
            % Crosstalk simulation
            if(flags.oct)
                
                %P_ct=1-exp(-(V_ob*gain/K_noise)^(1+Alpha_noise));
                
                Noct_seeds=poissrnd(lambda_prompt_ct_operative*mu);
                
                Ioct=[];
                Joct=[];
                I_neighboring_cell=[i+1,i,i-1,i];
                J_neighboring_cell=[j,j-1,j,j+1];
    
                for n=1:Noct_seeds   
                    Celda=ceil(rand(1)*4);
                    Ioct=[Ioct,I_neighboring_cell(Celda)];
                    Joct=[Joct,J_neighboring_cell(Celda)];
        
                end
                
                for mp = 1:Noct_seeds,
                   
                    ip=Ioct(mp);
                    jp=Joct(mp);
                    

                    time_oct = mcell(i, j,ind_t).t(k) + ( mp/10)*oct_delay;
                    if(ip > 0 && ip <= sipm_pix_xno && jp > 0 && jp <= sipm_pix_yno && time_oct < tmeas)
                        ind_ttmp = ind_t;
                        
                        for tmp_ind_t = ind_t+1:num_time_samples,
                            if( time_oct < Tsamp(tmp_ind_t))
                                ind_ttmp = tmp_ind_t-1;
                                break;
                            end
                        end
                        if(time_oct > Tsamp(num_time_samples) && time_oct < tmeas)
                            ind_ttmp = num_time_samples;
                        end
                        
                        mcell(ip,jp,ind_ttmp ).t = [mcell(ip,jp,ind_ttmp).t , mcell(i, j,ind_t).t(k) + ( mp/10)*oct_delay];
                        mcell(ip,jp,ind_ttmp ).count = mcell(ip,jp,ind_ttmp ).count + 1;
                        mcell(ip,jp,ind_ttmp).gain = [mcell(ip,jp,ind_ttmp).gain 0];
                        [mcell(ip,jp,ind_ttmp).t, order] = sort(mcell(ip,jp,ind_ttmp).t,'ascend');
                        %                                 mcell(ip,
                        %                                 jp,ind_t).flag = mcell(ip,jp,ind_t).flag(order);
                        %                             mcell(ip,jp,ind_ttmp ).flag = [mcell(ip,jp,ind_ttmp ).flag 4];
                        %                                 kp = mcell(ip,jp,ind_t).count;
                        
                        phot_pot.t = [phot_pot.t time_oct];
                        phot_pot.flag = [phot_pot.flag 4];
                        phot_pot.i = [phot_pot.i ip];
                        phot_pot.j = [phot_pot.j jp];
                        phot_pot.count = phot_pot.count + 1;
                        phot_pot.gain = [phot_pot.gain 0];
                        phot_pot_samp_count(ind_ttmp) = phot_pot_samp_count(ind_ttmp) + 1;              
                        
                        if mcell_hit_first_gains(i,j).first_photon==0 && (phot_pot.flag(phot_glob_count)==1 || mcell_hit_first_gains(i,j).noise_count<=length(mcell_hit_first_gains(i,j).first_noise_time))
                            if phot_pot.flag(phot_glob_count)~=1
                                
                                i_origin=  mcell_hit_first_gains(i,j).first_noise_origin_i(mcell_hit_first_gains(i,j).noise_count);
                                j_origin=  mcell_hit_first_gains(i,j).first_noise_origin_j(mcell_hit_first_gains(i,j).noise_count);

                                mcell_hit_first_gains(ip,jp).first_noise_origin_i=[mcell_hit_first_gains(ip,jp).first_noise_origin_i i_origin];
                                mcell_hit_first_gains(ip,jp).first_noise_origin_j=[mcell_hit_first_gains(ip,jp).first_noise_origin_j j_origin];
                                mcell_hit_first_gains(ip,jp).first_noise_time=[mcell_hit_first_gains(ip,jp).first_noise_time time_oct];
                            else
                                mcell_hit_first_gains(ip,jp).first_noise_origin_i=[mcell_hit_first_gains(ip,jp).first_noise_origin_i i];
                                mcell_hit_first_gains(ip,jp).first_noise_origin_j=[mcell_hit_first_gains(ip,jp).first_noise_origin_j j];
                                mcell_hit_first_gains(ip,jp).first_noise_time=[mcell_hit_first_gains(ip,jp).first_noise_time time_oct];

                            end

                            [mcell_hit_first_gains(ip,jp).first_noise_time,order] = sort(mcell_hit_first_gains(ip,jp).first_noise_time, 'ascend');
                            mcell_hit_first_gains(ip,jp).first_noise_origin_i = mcell_hit_first_gains(ip,jp).first_noise_origin_i(order);
                            mcell_hit_first_gains(ip,jp).first_noise_origin_j = mcell_hit_first_gains(ip,jp).first_noise_origin_j(order);
       

                            
                        end
                    end
                end
            end
            %
            % VM: Delayed Crosstalk simulation
            
            if(flags.doct)
                
                %prob_ct=1-exp(-(V_ob*gain/K_noise)^(1+Alpha_noise));
                
                Noct = poissrnd(lambda_delayed_ct_operative*gain);
                Ioct=[];
                Joct=[];
                I_neighboring_cell=[i+1,i,i-1,i];
                J_neighboring_cell=[j,j-1,j,j+1];
    
                for n=1:Noct   
                    Celda=ceil(rand(1)*4);
                    Ioct=[Ioct,I_neighboring_cell(Celda)];
                    Joct=[Joct,J_neighboring_cell(Celda)];
        
                end

                for mp = 1:Noct,
        
                    ip=Ioct(mp);
                    jp=Joct(mp);
        
                    Delay_dct=get_delayed_noise_time(1,P_dct,T_dct);
                    time_oct = mcell(i, j,ind_t).t(k) + Delay_dct;
                    if(ip > 0 && ip <= sipm_pix_xno && jp > 0 && jp <= sipm_pix_yno && time_oct < tmeas)
                        ind_ttmp = ind_t;
                        
                        for tmp_ind_t = ind_t+1:num_time_samples,
                            if( time_oct < Tsamp(tmp_ind_t))
                                ind_ttmp = tmp_ind_t-1;
                                break;
                            end
                        end
                        if(time_oct > Tsamp(num_time_samples) && time_oct < tmeas)
                            ind_ttmp = num_time_samples;
                        end
                        
                        mcell(ip,jp,ind_ttmp ).t = [mcell(ip,jp,ind_ttmp).t , mcell(i, j,ind_t).t(k) + Delay_dct];
                     	mcell(ip,jp,ind_ttmp ).count = mcell(ip,jp,ind_ttmp ).count + 1;
                        mcell(ip,jp,ind_ttmp).gain = [mcell(ip,jp,ind_ttmp).gain 0];
                        [mcell(ip,jp,ind_ttmp).t, order] = sort(mcell(ip,jp,ind_ttmp).t,'ascend');
                        %                                 mcell(ip,
                        %                                 jp,ind_t).flag = mcell(ip,jp,ind_t).flag(order);
                        %                             mcell(ip,jp,ind_ttmp ).flag = [mcell(ip,jp,ind_ttmp ).flag 4];
                        %                                 kp = mcell(ip,jp,ind_t).count;
                        
                        phot_pot.t = [phot_pot.t time_oct];
                        phot_pot.flag = [phot_pot.flag 5];
                        phot_pot.i = [phot_pot.i ip];
                        phot_pot.j = [phot_pot.j jp];
                        phot_pot.count = phot_pot.count + 1;
                        phot_pot.gain = [phot_pot.gain 0];
                        phot_pot_samp_count(ind_ttmp) = phot_pot_samp_count(ind_ttmp) + 1;              
                        
                         if mcell_hit_first_gains(i,j).first_photon==0 && (phot_pot.flag(phot_glob_count)==1 || mcell_hit_first_gains(i,j).noise_count<=length(mcell_hit_first_gains(i,j).first_noise_time))
                            if phot_pot.flag(phot_glob_count)~=1
                                
                                i_origin=  mcell_hit_first_gains(i,j).first_noise_origin_i(mcell_hit_first_gains(i,j).noise_count);
                                j_origin=  mcell_hit_first_gains(i,j).first_noise_origin_j(mcell_hit_first_gains(i,j).noise_count);

                                mcell_hit_first_gains(ip,jp).first_noise_origin_i=[mcell_hit_first_gains(ip,jp).first_noise_origin_i i_origin];
                                mcell_hit_first_gains(ip,jp).first_noise_origin_j=[mcell_hit_first_gains(ip,jp).first_noise_origin_j j_origin];
                                mcell_hit_first_gains(ip,jp).first_noise_time=[mcell_hit_first_gains(ip,jp).first_noise_time time_oct];
                            else
                                mcell_hit_first_gains(ip,jp).first_noise_origin_i=[mcell_hit_first_gains(ip,jp).first_noise_origin_i i];
                                mcell_hit_first_gains(ip,jp).first_noise_origin_j=[mcell_hit_first_gains(ip,jp).first_noise_origin_j j];
                                mcell_hit_first_gains(ip,jp).first_noise_time=[mcell_hit_first_gains(ip,jp).first_noise_time time_oct];

                            end
                             
                            
                            
                            [mcell_hit_first_gains(ip,jp).first_noise_time,order] = sort(mcell_hit_first_gains(ip,jp).first_noise_time, 'ascend');
                            mcell_hit_first_gains(ip,jp).first_noise_origin_i = mcell_hit_first_gains(ip,jp).first_noise_origin_i(order);
                            mcell_hit_first_gains(ip,jp).first_noise_origin_j = mcell_hit_first_gains(ip,jp).first_noise_origin_j(order);
       

                            
                        end
           
                    end
                end
            end

            
            if mcell_hit_first_gains(i,j).noise_count<=length(mcell_hit_first_gains(i,j).first_noise_time)  && mcell_hit_first_gains(i,j).first_noise_time( mcell_hit_first_gains(i,j).noise_count )==mcell(i,j,ind_t).t(k)
                i_origin=  mcell_hit_first_gains(i,j).first_noise_origin_i(mcell_hit_first_gains(i,j).noise_count);
                j_origin=  mcell_hit_first_gains(i,j).first_noise_origin_j(mcell_hit_first_gains(i,j).noise_count);

                mcell_hit_first_gains(i_origin,j_origin).first_gain= mcell_hit_first_gains(i_origin,j_origin).first_gain+gain;
                mcell_hit_first_gains(i,j).noise_count=mcell_hit_first_gains(i,j).noise_count+1;
            end

            if mcell_hit_first_gains(i,j).first_photon==0 && phot_pot.flag(phot_glob_count)==1
                mcell_hit_first_gains(i,j).first_photon=1;
                mcell_hit_first_gains(i,j).first_gain=mcell_hit_first_gains(i,j).first_gain+gain;
            end

        end
        %
        
        %%%
        
        % VM: sorting events
        [phot_pot.t,order] = sort(phot_pot.t, 'ascend');
        phot_pot.flag = phot_pot.flag(order);
        phot_pot.i = phot_pot.i(order);
        phot_pot.j = phot_pot.j(order);
        phot_pot.gain = phot_pot.gain(order);
        
        
        % VM: checking if the is another event on current time sample
        if(phot_count == phot_pot_samp_count(ind_t)) % All the photons have been processed. Using a while loop due to afterpulsing and crosstalk, which can
            %create new photons within the same time interval as well
            photon_present = 0;
            phot_pot.analysed_count = phot_pot.analysed_count + phot_pot_samp_count(ind_t);
        end
        %
        
    end
    %%%
    
    % VM:Calculus of the Current intensity through the Rs load
    if flags.simplified_readout == 0
        I_Rs = 0;
        if(flags.output_voltage)
            if(ind_t < num_time_samples)
                for k = 1:phot_pot.analysed_count,
                    tp = t_samp + delt(ind_t) - phot_pot.t(k);
                    I_Rs = I_Rs + G*(charge_elec) * phot_pot.gain(k) * ...
                    ( a1*exp(-tp/taucd1) + a2* exp(-tp/taucd2) +  a3 *exp(-tp/taumd) + a4 * exp(-tp/taumr));
                end
            end
        end
    end
    %     charge_time_int(ind_t) = I_Rs*delt(ind_t)/(G*charge_elec);
    %
    ind_t = ind_t+1;
end

%%%

%%% VM:Adding the created charge

%ph_fired = sum( phot_pot.gain)*G*charge_elec; % Final charge in coulombs
ph_fired = sum( phot_pot.gain); % Final charge in gain units

Q_first=[];
for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        Q_first=[Q_first, mcell_hit_first_gains(i,j).first_gain];
    end
end
phot_pot.first_gain_mean=Q_first;


