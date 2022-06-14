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
dtscale = 10; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);


% Set flow rate and cost to switch
Qave = V(end)/t(end);
Q = Qave*5; % Arbitraily set flow rate - to be changed later
Qp = [0 Q];
Cost2Switch = [0;1e-4];
%%
% Lets say that at the end of time, the accumulator is at zero net flow.
Vp = V(end);
%for k =1:length(t)
k = 1;
Flows = unique(Vp(1:k,k),'stable')-Qp*DPdt;
Vp(1:2*k,k+1) = Flows(:)

Cost(1:2*k,k) = abs(Vp(1:2*k,k+1)-V(end-dtscale*k)) + Cost2Switch

k = 2;
Flows = unique(Vp(1:2*(k-1),k),'stable')-Qp*DPdt;
Vp(1:2*k,k+1) = Flows(:)

Cost2Switch = Cost2Switch + Cost2Switch.';
Cost2Switch = Cost2Switch(:);

Cost(1:2*k,k) = abs(Vp(1:2*k,k+1)-V(end-dtscale*k)) + Cost2Switch

k = 3;
Flows = unique(Vp(1:2*(k-1),k),'stable')-Qp*DPdt;
Vp(1:2*k,k+1) = Flows(:)

Cost2Switch = Cost2Switch + Cost2Switch.';
Cost2Switch = Cost2Switch(:);

Cost(1:2*k,k) = abs(Vp(1:2*k,k+1)-V(end-dtscale*k)) + Cost2Switch

%%
for k = 1:3%length(DPt)
    Flows = Vp(1:k,k)-Qp*DPdt;
    Vp(1:2*k,k+1) = Flows(:);

    Cost2Switch = Cost2Switch + Cost2Switch.';
    Cost2Switch = Cost2Switch(:);

    Cost(1:2*k,k) = abs(Vp(1:2*k,k+1)-V(end-dtscale*k)) + Cost2Switch;

end



