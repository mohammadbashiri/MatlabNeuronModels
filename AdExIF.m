%% Addaptive Exponential I&F Model

close all

% initializing simulation param
tend = 1000;
fs   = 100;
dt   = 1/fs;
t    = 0:dt:tend-dt;
N    = numel(t);

% Stimulation current
I_stim = 1000;
I = [ones(1, 0.1*N)*0 ones(1, 0.5*N)*I_stim ones(1, 0.4*N)*0]; %ones(1,N)*600;

% initializations
C     = 281;    % pF
gL    = 30;     % nS
EL    = -70.6;  % mV
VT    = -50.4;  % mV
dT    = 2;      % mV
tauw  = 144;    % ms
a     = 4;      % nS
b     = 50; %0.0805; % nA
vpeak = 20;     % mV

v = ones(1,N) * EL;
w = zeros(1,N);

for i=1:N-1
    
    f    = -gL * (v(i) - EL) + gL * dT * exp((v(i) - VT)/dT);
    dvdt = ((f - w(i) + I(i))/C) * dt;
    v(i+1) = v(i) + dvdt;
    
    dwdt = ((a * (v(i) - EL) - w(i)) / tauw) * dt;
    w(i+1) = w(i) + dwdt;
    
    if v(i+1) > vpeak
        v(i+1) = EL;
        w(i+1) = w(i+1) + b;
    end
end


figure(1);
subplot(5,1,[1, 2]); plot(t, I); ylim([min(I)-10,max(I)+10]); grid;
legend('Current Density', 'Location','northeast');
ylabel({'$I(pA)$'},'Interpreter','latex');

subplot(5,1,[3, 4, 5]); plot(t, v); ylim([-100, 60]); grid;
legend('Action Potential','Location','northeast')
xlabel({'$Time (ms)$'},'Interpreter','latex');
ylabel({'$V_m (mV)$'},'Interpreter','latex');

suptitle({'Hodgkin Huxley Model', '(Giant Squid Axon)'});