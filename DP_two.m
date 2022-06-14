%% Deep Programming
% Fill / empty the rails so that the accumulators sizes are minimized.

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')
%figure(1), plot(t,cumsum(QR_1),t,cumsum(QR_2),t,cumsum(QR_3))
%legend('Rail 1','Rail 2','Rail 3'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')

% Posative flow is flow leaving the accumulator

% We are only pumping to rail 3
V = cumsum(QR_3);

% Start at t(end)
dt = t(2)-t(1);
dtscale = 1000; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);


% Set flow rate and cost to switch
Qave = V(end)/t(end);
Q = 2*Qave; % Arbitraily set flow rate - to be changed later
Qp = [0;Q];
Cost2Switch = 1e-9;
%
% Lets say that at the end of time, the accumulator is at zero net flow.
Vp = V(end);
Cost = 0;
for k =1:length(DPt)-1
    % flows at next stage
    ind = 1;
    for i = 1:size(Vp(:,k),1)
        Vp(ind:ind+1,k+1) = Vp(i,k) - Qp*DPdt;
        ind = ind+2;
    end
    
    % costs to get to that stage
    for i = 1:size(Vp(:,k+1),1)
        %cost due to error
        Error_cost =  abs( Vp(i,k+1) - V(find(t== DPt(length(DPt)-k))) );
        
        % Cost due to switching
        if ((rem(i,2) ==0) + (rem(round(i/2),2) ==0)) == 1
            Switching_cost = Cost2Switch;
        else
            Switching_cost = 0;
        end
        
        Cost(i,k+1) = Cost(round(i/2),k) + Error_cost + Switching_cost;
    end
    

    % If 2 options are identical, ie, they have the same cost and flow,
    % then only keep one
    [UniqueVp,ind_Vp] = unique(Vp(:,k+1));
    [UniqueCost,ind_Cost] = unique(Cost(:,k+1));
end
[mincost,ind_mincost] = min(Cost(:,end));

%Get back the optimal trajectory
ind_sequence = ind_mincost;
for i = 1:length(DPt)-1
    ind_sequence = round(ind_sequence/2);
    IsPumping(length(DPt)+1-i) = rem(ind_sequence,2) == 0;
end

ind = ind_mincost;
for i = 1:length(DPt)
    PumpFlow(i) = Vp(ind,end-i+1);
    ind = round(ind/2);
end
figure(3)
plot(DPt,PumpFlow,t,V), legend('Flow Provided by Pump','Flow Required','location','NorthWest')

%check that the cost checks out
TotalCostCheck = abs(PumpFlow(1)-V(find(t==DPt(1))) ); % Cost at beginning of time - to initize the variable
for i = 2:length(DPt)
    TotalCostCheck = TotalCostCheck + abs(PumpFlow(i)-V(find(t==DPt(i))) ) + Cost2Switch*abs(IsPumping(i)-IsPumping(i-1));
end
[TotalCostCheck; mincost]
DPt(end)
