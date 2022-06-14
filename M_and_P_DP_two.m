clear
load JCB5T_C0P_3CPR_Flows.mat
dt = t(2)-t(1);
inst_flow = [QR_1 QR_2 QR_3];
cumulative_flow = cumsum(inst_flow)*dt;
%figure(), plot(t,cumulative_flow), legend('Rail 1','Rail 2','Rail 3'), ylabel('Volume (m^3'), xlabel('Time (s)')

dtscale = 100;
DPdt = dt*dtscale;
DPtime = 0:DPdt:t(end);

Vpump = cumulative_flow(:,3); % Flow used by the actuators
Vmotor = cumulative_flow(:,1); % Flow used by the actuators


Qmin = max(Vpump(end) ,-Vmotor(end))/t(end);
Qscale = 20;
Q = Qscale*Qmin; % fixed flow rate of main pump/motor

Vpump_supply_options = 0;
Vmotor_supply_options = 0;
for i = 1: Qmin*t(end)/Q/DPdt + 5
    Vpump_supply_options = [Vpump_supply_options; Q*DPdt*i]; % Options for the pumping rail
    Vmotor_supply_options = [Vmotor_supply_options; -Q*DPdt*i]; % Options for the motoring rail
end

[Voptions1,Voptions2] = meshgrid(Vpump_supply_options, Vmotor_supply_options);
Voptionspump = Voptions1(:);
Voptionsmotor = Voptions2(:);

close all
figure()
plot(t,Vpump,t,Vmotor), hold on
for ii = 1:length(DPtime)
    %PercentComplete = ii/length(DPtime)*100
    % The four options are:
    
    % 1) motor the motoring rail by connecting it to tank (pump rail is not
    % affected)
    MT_pump = Voptionspump;
    for i = 1:length(Voptionspump)
        if Voptionsmotor(i) == min(Voptionsmotor)
            MT_motor(i) = inf;
        else
            MT_motor(i) = Voptionsmotor(i) - Q*DPdt;
        end
    end
    
    % 2) motor the motoring rail by pumping to the pumping rail
    for i = 1:length(Voptionspump)
        if Voptionsmotor(i) == min(Voptionsmotor)
            MP_motor(i) = inf;
            MP_pump(i) = inf;
        elseif Voptionspump(i) == max(Voptionspump)
            MP_pump(i) = inf;
            MP_motor(i) = inf;
        else
            MP_motor(i) = Voptionsmotor(i) - Q*DPdt;
            MP_pump(i) = Voptionspump(i) + Q*DPdt;
        end
    end
    
    
    % 3) pump the pumping rail by connecting it to tank (motoring rail is
    % unaffected)
    TP_motor = Voptionsmotor;
    for i = 1:length(Voptionspump)
        if Voptionspump(i) == max(Voptionspump)
            TP_pump(i) = inf;
        else
            TP_pump(i) = Voptionspump(i) + Q*DPdt;
        end
    end
    
    % 4) pump tank to tank
    TT_motor = Voptionsmotor;
    TT_pump = Voptionspump;
    
    
    %% Which option was best at each point?
    [~,t_ind_forward] = min(abs(t - t(end)+(ii-1)*DPdt));
    [~,t_ind_backward] = min(abs(t - t(end)+ii*DPdt));
    
    Error1 = (abs(MT_motor'-Vmotor(t_ind_forward))+abs(MT_pump-Vpump(t_ind_forward)));
    Error2 = (abs(MP_motor'-Vmotor(t_ind_forward))+abs(MP_pump'-Vpump(t_ind_forward)));
    Error3 = (abs(TP_motor-Vmotor(t_ind_forward))+abs(TP_pump'-Vpump(t_ind_forward)));
    Error4 = (abs(TT_motor-Vmotor(t_ind_forward))+abs(TT_pump-Vpump(t_ind_forward)));
    [cost,thischoice] = min([Error1  Error2  Error3  Error4]');
    choices(:,ii) = thischoice';
    
    for j = 1:length(choices)
        if choices(j,ii) = 1
    end
    
    clear MT_motor MT_pump MP_motor MP_pump TP_motor TP_pump TT_motor TT_pump
end


%% Get optimal trajectory back
choice(1