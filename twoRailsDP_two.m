%% Dynamic Programming THIS IS JAKE's algorithm
% Fill / empty the rails so that the accumulators sizes are minimized.
return
% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.

% This code takes a drive cycle, and modifies it so that it has two rails
% that need main pump flow

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')
% figure(1), plot(t,cumsum(QR_1),t,cumsum(QR_2),t,cumsum(QR_3))
% legend('Rail 1','Rail 2','Rail 3'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')

% Posative flow is flow leaving the accumulator

% We are pumping to rail 3 and inverting rail 1 and pumping to it
V1 = cumsum(QR_3);
V2 = cumsum(-QR_1/2);
V = [V1 V2];
% figure(1), plot(t,V1,t,V2)
% legend('Rail 1','Rail 2'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')


% Start at t(end)
dt = t(2)-t(1);
dtscale = 100; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);

% Set flow rate and cost to switch
Qave = max(V(end,:)/t(end)) ;
Q = 2*Qave; % Arbitraily set flow rate - to be changed later
Cost2Switch = 0*1e-2;

% nn is the number of time steps with flow it takes to get the the required
% volume
nn = max(ceil(V1(end)/Q/DPdt) +2,ceil(V1(end)/Q/DPdt) +2);


% Make V_MP - which is a matrix of possible values of V_p at each time
for i = 0:nn-1
V_MP1(i+1) = i*Q*DPdt;
end
V_MP = repmat(V_MP1,1,size(V,2));


% Build the cost
J = NaN(nn^size(V,2),length(DPt));
ind = J;

% Make indexers
a = fliplr(cumsum(ones(1,length(V_MP1)))); b = 2*a(1):-1:a(1)+1; c = 3*a(1):-1:2*a(1)+1;
if size(V,2) == 3 
    [x,y,z] = meshgrid(a,b,c); X = x(:); Y = y(:); Z = z(:);
    J(:,end) = abs(V(end,1)-V_MP(X)) + abs(V(end,2)-V_MP(Y)) + abs(V(end,end)-V_MP(Z)) ;
    disp('Size of V not accounted for in this code')
    return
elseif size(V,2) == 2
    [x,y] = meshgrid(a,b); X = x(:); Y = y(:);
    J(:,end) = abs(V(end,1)-V_MP(X)) + abs(V(end,2)-V_MP(Y));
elseif size(V,2) ==1
    J(:,end) = abs(V(end,1)-V_MP);
else
    disp('Size of V not accounted for in this code')
    return
end
%



%
%PercentDone = 0
for k = 1:length(DPt)-1
    [~,t_ind] = min(abs(t-(t(end)-k*DPdt)));
    J(1,end-k) = abs(V1(t_ind)-V_MP(1))+ abs(V1(t_ind)-V_MP(1)) + J(1,end-k+1);
    ind(1,end-k) = 1;
    for j = 1:nn^size(V,2)-1
        if X(j+1) == max(X) % The case where rail 1 should not be filled anymore
            [J(j+1,end-k),ind(j+1,end-k)] = min([J(j+1,end-k+1)+Cost2Switch*(ind(j+1,end-k+1)==2|ind(j+1,end-k+1)==3),inf,J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2)]);
            J(j+1,end-k) = J(j+1,end-k) + abs(V(t_ind,1)-V_MP(X(j+1))) + abs(V(t_ind,2)-V_MP(Y(j+1)));
            
        elseif Y(j+1) == max (Y) % The case where rail 2 should not be filled anymore
            [J(j+1,end-k),ind(j+1,end-k)] = min([J(j+1,end-k+1)+Cost2Switch*(ind(j+1,end-k+1)==2|ind(j+1,end-k+1)==3),J(j+1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j+1- length(V_MP1),end-k+1)==1 | ind(j+1- length(V_MP1),end-k+1)==3),inf]);
            J(j+1,end-k) = J(j+1,end-k) + abs(V(t_ind,1)-V_MP(X(j+1))) + abs(V(t_ind,2)-V_MP(Y(j+1)));            
            
        else
            [J(j+1,end-k),ind(j+1,end-k)] = min([J(j+1,end-k+1)+Cost2Switch*(ind(j+1,end-k+1)==2 | ind(j+1,end-k+1)==3),J(j+1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j+1- length(V_MP1),end-k+1)==1 | ind(j+1- length(V_MP1),end-k+1)==3),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1 | ind(j,end-k+1)==2)]);
            J(j+1,end-k) = J(j+1,end-k) + abs(V(t_ind,1)-V_MP(X(j+1))) + abs(V(t_ind,2)-V_MP(Y(j+1))); 
            
        end

    end
    
    
    
    PercentDone = 100*k/length(DPt);
    if rem(PercentDone,10) < 100/length(DPt)
        PercentDone = (round(PercentDone,-1))
    end
end
PercentDone = 100


% Get optimal solution back out
% We have to start at zero flow
MinCost= min(J(end,1));

current_ind = size(J,1);
for ii = 1:length(DPt)-1
    if ind(current_ind,ii) == 1
        IsPumping(ii) = 0;
        current_ind = current_ind;
    elseif ind(current_ind,ii) == 2
        IsPumping(ii) = 1;
        current_ind = current_ind - length(V_MP1);
    elseif ind(current_ind,ii) == 3
        IsPumping(ii) = 2;
        current_ind = current_ind - 1;
    end
end

% Find what the volume looks like at each time
V_p1 = 0;
V_p2 = 0;
for jj = 1:length(DPt)-1
    V_p1(jj+1) = V_p1(jj) + (IsPumping(jj)==1)*Q*DPdt;
    V_p2(jj+1) = V_p2(jj) + (IsPumping(jj)==2)*Q*DPdt;
end



% Does the cost match?
CostCheck = abs(V(1,1)-V_p1(1)) + abs(V(1,2)-V_p2(1));
IsPumping(end+1) = IsPumping(end);
for kk  = 2:length(DPt)
    [~,t_ind] = min(abs( DPt(kk) - t ));
    CostCheck = CostCheck +  abs(V(t_ind,1)-V_p1(kk)) + abs(V(t_ind,2)-V_p2(kk));
    CostCheck = CostCheck + (IsPumping(kk-1)~=IsPumping(kk))*Cost2Switch;
end
[MinCost,CostCheck] %Those two numbers should be the same

% What is the required accumulator size?
for i = 1:length(t)
    [~,DPt_ind] = min(abs( DPt - t(i) ));
    Error1(i) = abs( V(i,1) - V_p1(DPt_ind));
    Error2(i) = abs( V(i,2) - V_p2(DPt_ind));
end
Accumulator_size1 = max(Error1); % m^3
Accumulator_size2 = max(Error2); % m^3
Accumulator_sizes = [Accumulator_size1 Accumulator_size2], disp('m^3')

figure(1)
plot(t,V,DPt,V_p1,DPt,V_p2), legend('Rail 1 Flow','Rail 2 Flow','Flow delivered to Rail 1 by Pump','Flow delivered to Rail 2 by Pump','Location','NorthWest')
title(['Accumulator Sizes = [',num2str(Accumulator_size1),' , ',num2str(Accumulator_size2) , ']',' m^3'])
ylabel('Flow m^3'), xlabel('Time(s)')
