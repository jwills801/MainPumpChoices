%% Dynamic Programming THIS IS JAKE's algorithm
% Fill / empty the rails so that the accumulators sizes are minimized.

% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')


% Start at t(end)
dt = t(2)-t(1);
dtscale = 10; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);

% figure(1), plot(t,cumsum(QR_1)*dt,t,cumsum(QR_2)*dt,t,cumsum(QR_3)*dt)
% legend('Rail 1','Rail 2','Rail 3'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')

% Posative flow is flow leaving the accumulator

% We are only pumping to rail 3
V = cumsum(QR_3)*dt;

% Set flow rate and cost to switch
Qave = V(end)/t(end);
N = 500;
Qvals = linspace(Qave,20*Qave,N);
Cost2Switch = 0*1e-1;

PercentDone = 0
for jjj = 1:length(Qvals)
    Q = Qvals(jjj);
    
    % nn is the number of time steps with flow it takes to get the the required
    % volume
    nn = ceil(V(end)/Q/DPdt) +2;
    
    % Make V_MP - which is a matrix of possible values of V_p at each time
    V_MP = 0;
    for i = 0:nn-1
        V_MP(i+1) = i*Q*DPdt;
    end
    V_MP = repmat(fliplr(V_MP)',1,length(DPt));
    
    %
    % Build the cost
    J = NaN(nn,length(DPt));
    ind = J;
    J(:,end) = abs(V(end)-V_MP(:,end));
    
    
    for k = 1:length(DPt)-1
        t_ind = find(abs(t-(t(end)-k*DPdt)) < min(dt/2,DPdt/2),1);
        J(1,end-k) = abs(V(t_ind)-V_MP(1,end-k+1)) + J(1,end-k+1);
        ind(1,end-k) = 1;
        for j = 1:nn-1
            [J(j+1,end-k),ind(j+1,end-k)] = min([J(j+1,end-k+1)+Cost2Switch*(ind(j+1,end-k+1)==2),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1)]);
            J(j+1,end-k) = J(j+1,end-k) + abs(V(t_ind)-V_MP(j+1,end-k));
        end
        
        
        

    end
    PercentDone = 100*jjj/length(Qvals);
    if rem(PercentDone,10) < 100/length(Qvals)
        PercentDone = (round(PercentDone,-1))
    end
    
    
    % Get optimal solution back out
    [MinCost,min_ind] = min(J(:,1));
    
    current_ind = min_ind;
    for ii = 1:length(DPt)
        if ind(current_ind,ii) == 1
            IsPumping(ii) = 0;
            current_ind = current_ind;
        elseif ind(current_ind,ii) == 2
            IsPumping(ii) = 1;
            current_ind = current_ind - 1;
        end
    end
    
    % Find what the volume looks like at each time
    V_p = 0;
    for jj = 1:length(DPt)-1
        V_p(jj+1) = V_p(jj) + IsPumping(jj)*Q*DPdt;
    end
    
    % Does the cost match?
    CostCheck = abs(V(1,1)-V_p(1));
    IsPumping(end+1) = IsPumping(end);
    for kk  = 2:length(DPt)
        [~,t_ind] = min(abs( DPt(kk) - t ));
        CostCheck = CostCheck +  abs(V(t_ind)-V_p(kk));
        CostCheck = CostCheck + (IsPumping(kk-1)~=IsPumping(kk))*Cost2Switch;
    end
    [MinCost,CostCheck]; %Those two numbers should be the same
    
    % What is the required accumulator size?
    for i = 1:length(t)
        [~,DPt_ind] = min(abs( DPt - t(i) ));
        Error(i) = abs( V(i) - V_p(DPt_ind));
    end
    Accumulator_size(jjj) = max(Error); % m^3
    %disp('m^3')
    
    
end


figure(2)
plot(Qvals./Qave,Accumulator_size), ylabel('Minimum Accumulator Size (m^3)'), xlabel('Flow Rate (multiples of average flow rate)'), title(['Minimum accumulator size is ',num2str(min(Accumulator_size)),' m^3'])
grid



%% BEST TIME STEP
%% Dynamic Programming THIS IS JAKE's algorithm
% Fill / empty the rails so that the accumulators sizes are minimized.

% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')

dtscalevals = 0:5:100; dtscalevals(1) = 1;

PercentDone = 0
for jjj = 1:length(dtscalevals)
    % Start at t(end)
    dt = t(2)-t(1);
    dtscale = dtscalevals(jjj); % Step through DP at a different time step than the one given by Drive Cycle
    DPdt = dtscale*dt;
    DPt = 0:DPdt:t(end);
    
    % figure(1), plot(t,cumsum(QR_1)*dt,t,cumsum(QR_2)*dt,t,cumsum(QR_3)*dt)
    % legend('Rail 1','Rail 2','Rail 3'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')
    
    % Posative flow is flow leaving the accumulator
    
    % We are only pumping to rail 3
    V = cumsum(QR_3)*dt;
    
    % set cost to switch
    Cost2Switch = 0*1e-1;
    Qave = V(end)/t(end);
    Q = 5.5*Qave;
    
    % nn is the number of time steps with flow it takes to get the the required
    % volume
    nn = ceil(V(end)/Q/DPdt) +2;
    
    % Make V_MP - which is a matrix of possible values of V_p at each time
    V_MP = 0;
    for i = 0:nn-1
        V_MP(i+1) = i*Q*DPdt;
    end
    V_MP = repmat(fliplr(V_MP)',1,length(DPt));
    
    %
    % Build the cost
    J = NaN(nn,length(DPt));
    ind = J;
    J(:,end) = abs(V(end)-V_MP(:,end));
    
    
    for k = 1:length(DPt)-1
        t_ind = find(abs(t-(t(end)-k*DPdt)) < min(dt/2,DPdt/2),1);
        J(1,end-k) = abs(V(t_ind)-V_MP(1,end-k+1)) + J(1,end-k+1);
        ind(1,end-k) = 1;
        for j = 1:nn-1
            [J(j+1,end-k),ind(j+1,end-k)] = min([J(j+1,end-k+1)+Cost2Switch*(ind(j+1,end-k+1)==2),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1)]);
            J(j+1,end-k) = J(j+1,end-k) + abs(V(t_ind)-V_MP(j+1,end-k));
        end
        
        
        
        
    end
    PercentDone = 100*jjj/length(dtscalevals);
    if rem(PercentDone,10) < 100/length(dtscalevals)
        PercentDone = (round(PercentDone,-1))
    end
    
    
    % Get optimal solution back out
    [MinCost,min_ind] = min(J(:,1));
    
    current_ind = min_ind;
    for ii = 1:length(DPt)
        if ind(current_ind,ii) == 1
            IsPumping(ii) = 0;
            current_ind = current_ind;
        elseif ind(current_ind,ii) == 2
            IsPumping(ii) = 1;
            current_ind = current_ind - 1;
        end
    end
    
    % Find what the volume looks like at each time
    V_p = 0;
    for jj = 1:length(DPt)-1
        V_p(jj+1) = V_p(jj) + IsPumping(jj)*Q*DPdt;
    end
    
    % Does the cost match?
    CostCheck = abs(V(1,1)-V_p(1));
    IsPumping(end+1) = IsPumping(end);
    for kk  = 2:length(DPt)
        [~,t_ind] = min(abs( DPt(kk) - t ));
        CostCheck = CostCheck +  abs(V(t_ind)-V_p(kk));
        CostCheck = CostCheck + (IsPumping(kk-1)~=IsPumping(kk))*Cost2Switch;
    end
    [MinCost,CostCheck]; %Those two numbers should be the same
    
    % What is the required accumulator size?
    for i = 1:length(t)
        [~,DPt_ind] = min(abs( DPt - t(i) ));
        Error(i) = abs( V(i) - V_p(DPt_ind));
    end
    Accumulator_size(jjj) = max(Error); % m^3
    %disp('m^3')
end

figure
plot(dtscalevals*dt,Accumulator_size), ylabel('Minimum Accumulator Size (m^3)'), xlabel('Time step'), title(['Minimum accumulator size is ',num2str(min(Accumulator_size)),' m^3'])
grid
















