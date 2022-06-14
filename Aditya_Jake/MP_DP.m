%% Goal: Find the optimal trajectory of main pump outlet decisions.

%% Assumptions

% 1. Assume there is a main pump for each rail.
% 2. Assume each MP is sized "SMP" times the minimum needed size. The
%    minimum size can be found by dividing net rail flow volume by the drive
%    cycle time and the assumed angular velocity.
    SMP = 1;
% 3. Drive cycle is ...
%     load RF_LC_5CPRs_CMOOP.mat
%     load RF_SC_5CPRs_CMOOP.mat
    load JCB5T_C0P_3CPR_Flows.mat


%% Cost Function

% The cost function will be the largest absolute accumulator volume change.
% The accumular volume change can be found by integrating the difference in
% HECM and MP flow rates.

%% Constraints

% 1. Net accumular volume at beginning and end of drive cycle must be
%    zero.
% 2. Note there is a maximum amount of feasible flow each time step and
%    that the main pump is either pumping fluid into the rail, or it is
%    not. Q_MP is in {0,dQ} where q is the main pump flow rate and dt
%    is the time step. dQ = q*dt.

%% Dynamic Programming

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE LATER
t = t(1:10000);
% Q1 = Q5(1:10000);
Q1 = QR_2(1:10000); %rail 2
Q2 = QR_3(1:10000); %rail 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The JCB Drive cycle has a 3 rails of which the middle one is constrained two have zero net flow
%So the pump does not need to operate on that rail and the accumulator will take care of the flows for it


% Time Step
dt = mean(diff(t));
% Find the change in accumulator volume due to the HECM ONLY.
V_H = cumtrapz(t,Q2); %change Q1 to Q2 if you want to analyse a different rail.


dV = SMP*V_H(end)/(t(end)-t(1))*dt;
% dV = SMP*V_H(end)/(t(end)-t(1))*t_space;

% V_H_2 = cumtrapz(t,Q2);
% dV_2 = SMP*V_H_2(end)/(t(end)-t(1))*dt;

% dV = max(dV_1,dV_2); %choosing which flow rate will be provided by pump


% Time Step Indicies
k_a = 1:length(t);

% States: Define the states as the net volume 
nn = round(V_H(end)/dV,0)+1;

V_MP = linspace(0,V_H(end),nn);
V_MP = repmat(V_MP,length(k_a),1);

% Cost Function Matrix Initialization
J = zeros(size(V_MP));
temp = zeros(1,size(V_MP,2));
temp( V_MP(end,:) ~= V_H(end) ) = inf;
J(end,:) = temp;

% Pump/Motor On/Off Indicator
ind = zeros(size(J));

% Volume Error Matrix Initialization
E = abs(V_MP-V_H);
% E_2 = abs(V_MP - V_H_2);


for i = 1:length(k_a)-1
    
    k = k_a(length(k_a)-i);
    
    for j = 1:size(J,2)
                
        if j<size(J,2)
            [J(k,j), ind(k,j)] = min([max([J(k+1,j), E(k,j)]), max([J(k+1,j+1), E(k,j)])]);
        else
            [J(k,j), ind(k,j)] = min([max([J(k+1,j), E(k,j)]), inf, E(k,j)]);
        end
    end
    CompletionPercentage = 100*i/(length(k_a)-1)
end

%% Analyze Choices

n = 1;

for i = 1:length(k_a)-1
    ind_forward(i) = ind(i,n);
    if ind_forward(i) ~= 1
        n = n+1;
    end
    
    i
    n
end

ind_forward(ind_forward==1) = 0; %no supply to rail 1
ind_forward(ind_forward~=0) = 1; 

% ind_forward_2(ind_forward==1) = 1;
% ind_forward_2(ind_forward~=0) = 0;

V_MP_sel = cumsum([ind_forward*dV 0])';
% V_MP_sel_2 = cumsum([ind_forward_2*dV 0])';




%% I added this
% Plot

figure(1)
plot(t,V_H, t,V_MP_sel)


figure(2)














