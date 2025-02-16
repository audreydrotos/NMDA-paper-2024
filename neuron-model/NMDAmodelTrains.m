function [t, v] = NMDAmodelTrains()

% WB model neuron with 2 synaptic current inputs
% WB model is a point neuron so units are in densities (per cm^2 for
% surface area)
% To change to a spherical neuron, multiply densities by neuron surface area
% Parameters are set in function modeleqs below

% set model params
g_syn1 = 0.15;
g_syn2 = 0*0.25;
g_l = 0.7;

% clear; clf;

t_final=1000; % duration of simulation in msec
tspan = [0 t_final];

% initial conditions for model variables
v0=-64;
h0=h_inf(v0);
n0=n_inf(v0);
s10=0;
s20=0;
ICs = [v0; h0; n0; s10; s20];

% calculate input times
inputs = [200 233.33 266.67 300 333.33];

% Call the solver ode15s to numerically simulate the model
% output: t = time stamp vector. 
% vars = matrix with column 1 = voltage v, 
% column 2 = h Na inactivation gate, 
% column 3 = n K-dr activation gate,
% column 4 = s synaptic current gate
options = odeset('MaxStep',1);
% [t,vars] = ode15s(@modeleqs, tspan, ICs, options);

% pass additional parameters using an anonymous function
[t, vars] = ode15s(@(t,vars) modeleqs(t, vars, inputs), tspan, ICs, options);

% vectors of model variables at each time step:
v = vars(:,1);
%h = vars(:,2);
%n = vars(:,3);
s1 = vars(:,4);
s2 = vars(:,5);

% determine spike times and interspike intervals
[peaks, indxs]=findpeaks(v,'MINPEAKHEIGHT',-20);
if ~isempty(indxs) 
    spiketimes=t(indxs);
else
    spiketimes = [];
end

% plot output
figure('Position', [0 0 900 682])
subplot(4,1,1);
plot(t,v,'Color', '#03045e','Linewidth',2);
set(gca,'Fontsize',16);
% xlabel('t [ms]','Fontsize',20); ylabel('v [mV]','Fontsize',20);
ylim([-65 -50])

subplot(4,1,2)
plot(t,s1*g_syn1,'Color','#0077b6','Linewidth',2)
set(gca,'Fontsize',16);
% xlabel('t [ms]','Fontsize',20); ylabel('ampa current gate','Fontsize',16);
ylim([0 0.2])

subplot(4,1,3)
plot(t,s2*g_syn2,'Color','#00b4d8','Linewidth',2)
set(gca,'Fontsize',16);
% xlabel('t [ms]','Fontsize',20); ylabel('nmda current gate','Fontsize',16);
ylim([0 0.2])

subplot(4,1,4)
plot(t,s1+s2,'Color','#90e0ef','Linewidth',2)
set(gca,'Fontsize',16);
ylim([0 2])
% xlabel('t [ms]','Fontsize',20); ylabel('summed current gate','Fontsize',16);

function dvarsdt = modeleqs(t,vars,inputs)
% WB model point neuron equations
% activation gating of Na current is instantaneous function of voltage

    % model variables
    v = vars(1);
    h = vars(2);
    n = vars(3);
    s1 = vars(4);
    s2 = vars(5);

    c=1; % membrane capacitance in microF/cm^2
    g_k=9;  % max conductance of K-dr current (mS/cm^2)
    g_na=35; % max conductance of Na current (mS/cm^2)
    % g_l=0.7; % max conductance of membrane leak current (mS/cm^2)
    v_k=-90; % reversal potential of K-dr current
    v_na=55; % reversal potential of Na current
    v_l=-65; % reversal potential of leak current
    
    i_ext=0;  % external applied current to neuron
    
    % pre-synaptic spikes
    % T = 1/(modFreq/1000);

    % T=50;       % period in msec of repetitive pre-synaptic spikes
    presyn_spike_width = 2; % (msec) if you choose a very slow synaptic rise time
    % constant, you may need to make presyn_spike_width longer to see
    % effects. this was initially 1 i changed it to 4. 
    % calculate times for input spikes
    
    % calculate the input spike right before the first spike, if there is
    % one
    prior_inputs = inputs(inputs < t);
    if ~isempty(prior_inputs)
        maximum = max(prior_inputs);
        diff = t - maximum;
        if diff <= presyn_spike_width
            q = 1; % during the spike
        else
            q = 0; % not during the spike
        end
    else
        q = 0; % before spike
    end
    % if mod(t,T) <= presyn_spike_width && t > 10.0
    %    q=1;
    % else
    %    q=0;
    % end

    % post-synaptic current 1 ampa current
    % g_syn1=0.400;   % max conductance  (mS/cm^2)
    tau_d1=32; tau_r1=1.3; % time constants for decay and rise of synaptic current (ms)
    
    % post-synaptic current 2 nmda currnet
    % g_syn2=0.100;   % max conductance  (mS/cm^2) % added this to function
    % tau_d2=115; tau_r2=70; % time constants for decay and rise of synaptic current (ms)
    tau_d2=65; tau_r2=3.5;

    dvdt = (g_k*n^4*(v_k-v) + g_na*m_inf(v)^3*h*(v_na-v) + ...
           g_l*(v_l-v) - g_syn1*s1*v - g_syn2*s2*v +i_ext)/c;
    dhdt = alpha_h(v)*(1-h)-beta_h(v)*h; 
    dndt = alpha_n(v)*(1-n)-beta_n(v)*n; 
    ds1dt = q*(1-s1)/tau_r1-s1/tau_d1;
    ds2dt = q*(1-s2)/tau_r2-s2/tau_d2;
    
    dvarsdt = [dvdt; dhdt; dndt; ds1dt; ds2dt];

end

end %% for function
