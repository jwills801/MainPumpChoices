% Fill / empty the rails so that the accumulators sizes are minimized.

% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.

% This code takes a drive cycle, and modifies it so that it has two rails
% that need main pump flow

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')
% Posative flow is flow leaving the accumulator




% Start at t(end)
dt = t(2)-t(1);
dtscale = 20; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);


V1 = cumsum(QR_3)*dt;
V2 = cumsum(QR_1)*dt;
V = [V1 V2];
figure, plot(t,V1,t,V2)
legend('Rail 1','Rail 2'), ylabel('Volume (m^3)'), xlabel('Time (s)')

% Set flow rate and cost to switch
Qave = max(abs(V(end,:)/t(end)));
%Q = 6.82*Qave; % if DPdt = .1 secs then best Flow rate = 6.82 *Qave
%Q = 5.65*Qave; % if DPdt = .2 secs then best Flow rate = 5.65 *Qave  => perc_thresh = .05
Q = 7.6826e-04 ; % Best if DPdt = .2 secs, K = 1.4 and perc_thresh = .1
%Q = 2*Qave;
Cost2Switch = 0*1e-2;

perc_threshold = .1;
K = 1.4;

% nn is the number of time steps with flow it takes to get the the required
% volume
nn = max(ceil(V1(end)/Q/DPdt) +2,ceil(V2(end)/-Q/DPdt) +2);


% Make V_MP - which is a matrix of possible values of V_p at each time
for i = 0:nn-1
    V_MP1(i+1) = i*Q*DPdt;
end
V_MP = [(V_MP1) (-V_MP1)];


% Build the cost
J = NaN(nn^size(V,2),length(DPt));
ind = J;
Accumulator_size1 = J ;
Accumulator_size2 = J ;

% Make indexers
a = fliplr(cumsum(ones(1,length(V_MP1)))); b = 2*a(1):-1:a(1)+1; c = 3*a(1):-1:2*a(1)+1;
if size(V,2) == 3 
    [x,y,z] = meshgrid(a,b,c); X = x(:); Y = y(:); Z = z(:);
    J(:,end) = abs(V(end,1)-V_MP(X)) + abs(V(end,2)-V_MP(Y)) + abs(V(end,end)-V_MP(Z)) ;
    disp('This size of V not accounted for in this code')
    return
elseif size(V,2) == 2
    [x,y] = meshgrid(a,b); X = x(:); Y = y(:);

    delV1 = -(V(end,1)-V_MP(X));
    delV2 = -(V(end,2)-V_MP(Y));
    Accumulator_size1(:,end) = ((1+sign(delV1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV1)*perc_threshold).^(1/K)-1) .* delV1;
    Accumulator_size2(:,end) = ((1+sign(delV2)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV2)*perc_threshold).^(1/K)-1) .* delV2;

    J(:,end) = Accumulator_size1(:,end) + Accumulator_size2(:,end) ;
elseif size(V,2) ==1
    J(:,end) = abs(V(end,1)-V_MP);
else
    disp('This size of V not accounted for in this code')
    return
end



%
%PercentDone = 0
for k = 1:length(DPt)-1
    [~,t_ind] = min(abs(t-(t(end)-k*DPdt)));
    
    delV1 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    delV2 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    A1 = max( ((1+sign(delV1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV1)*perc_threshold).^(1/K)-1) .* delV1 ) ;
    A2 = max( ((1+sign(delV2)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV2)*perc_threshold).^(1/K)-1) .* delV2 ) ;
    
    Accumulator_size1(1,end-k) = max([A1 Accumulator_size1(1,end-k+1)]) ;
    Accumulator_size2(1,end-k) = max([A2 Accumulator_size2(1,end-k+1)]) ;
    J(1,end-k) = Accumulator_size1(1,end-k+1) + Accumulator_size2(1,end-k+1);

    ind(1,end-k) = 4;
    for j = 2:nn^size(V,2)
        if V_MP(X(j)) == max(V_MP) % The case where rail 1 should not be filled anymore
            [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),inf,inf,J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            
        elseif V_MP(Y(j)) == min (V_MP) % The case where rail 2 should not be emptied anymore
            [J(j,end-k),ind(j,end-k)] = min([inf,inf,J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            %J(j,end-k) = max( [ J(j,end-k) , abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))) ] );            
            
        else
            [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),J(j-1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-1-length(V_MP1),end-k+1)==1|ind(j-1- length(V_MP1),end-k+1)==3|ind(j-1- length(V_MP1),end-k+1)==4),J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            %J(j,end-k) = max([J(j,end-k) , abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))) ] ); 
            %J(j,end-k) = J(j,end-k) + abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))); 
            
        end

        delV1 = -(V(t_ind:t_ind+dtscale-1,1)-V_MP(X(j)));
        delV2 = -(V(t_ind:t_ind+dtscale-1,2)-V_MP(Y(j)));
        A1 = max( ((1+sign(delV1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV1)*perc_threshold).^(1/K)-1) .* delV1 ) ;
        A2 = max( ((1+sign(delV2)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV2)*perc_threshold).^(1/K)-1) .* delV2 ) ;
        if ind(j,end-k) == 1
            Accumulator_size1(j,end-k) = max([A1 Accumulator_size1(j-1,end-k+1)]) ;
            Accumulator_size2(j,end-k) = max([A2 Accumulator_size2(j-1,end-k+1)]) ;
        elseif ind(j,end-k) == 2
            Accumulator_size1(j,end-k) = max([A1 Accumulator_size1(j-1-length(V_MP1),end-k+1)]) ;
            Accumulator_size2(j,end-k) = max([A2 Accumulator_size2(j-1-length(V_MP1),end-k+1)]) ;
        elseif ind(j,end-k) == 3
            Accumulator_size1(j,end-k) = max([A1 Accumulator_size1(j-length(V_MP1),end-k+1)]) ;
            Accumulator_size2(j,end-k) = max([A2 Accumulator_size2(j-length(V_MP1),end-k+1)]) ;
        else
            Accumulator_size1(j,end-k) = max([A1 Accumulator_size1(j,end-k+1)]) ;
            Accumulator_size2(j,end-k) = max([A2 Accumulator_size2(j,end-k+1)]) ;
        end

        J(j,end-k) = Accumulator_size1(j,end-k) + Accumulator_size2(j,end-k);


    end
    
    
    
    PercentDone = 100*k/length(DPt);
    if rem(PercentDone,10) < 100/length(DPt)
        PercentDone = (round(PercentDone,-1))
    end
end
PercentDone = 100



% Get optimal solution back out
% We have to start at zero flow
MinCost = min(J(end,1))
Accumulator_size1(end,1) + Accumulator_size2(end,1)

current_ind(1) = size(J,1);
for ii = 1:length(DPt)-1
    if ind(current_ind(ii),ii) == 1
        current_ind(ii+1) = current_ind(ii)-1;
    elseif ind(current_ind(ii),ii) == 2
        current_ind(ii+1) = current_ind(ii) - 1 - length(V_MP1);
    elseif ind(current_ind(ii),ii) == 3
        current_ind(ii+1) = current_ind(ii) - length(V_MP1);
    elseif ind(current_ind(ii),ii) == 4
        current_ind(ii+1) = current_ind(ii);
    end
end

% Find what the volume looks like at each time
V_p1 = 0;
V_p2 = 0;
for jj = 1:length(DPt)-1
    V_p1(jj+1) = V_p1(jj) + (ind(current_ind(jj),jj)==2|ind(current_ind(jj),jj)==3)*Q*DPdt;
    V_p2(jj+1) = V_p2(jj) - (ind(current_ind(jj),jj)==1|ind(current_ind(jj),jj)==2)*Q*DPdt;
end


% What is the required accumulator size?
for i = 1:length(t)
    [~,DPt_ind] = min(abs( DPt - t(i) + DPdt/2 - t(2)));
    Error1(i) = V(i,1) - V_p1(DPt_ind);
    Error2(i) = V(i,2) - V_p2(DPt_ind);
end
DelV1 = -Error1 ;
DelV2 = -Error2 ;
figure, plot(t,DelV1), ylabel('$\Delta V (m^3)$','interpreter','latex'), xlabel('Time (s)','interpreter','latex')
figure, plot(t,DelV2), ylabel('$\Delta V (m^3)$','interpreter','latex'), xlabel('Time (s)','interpreter','latex')


delV_max1 = max(DelV1) ;
delV_min1 = min(DelV1) ;
delV_max2 = max(DelV2) ;
delV_min2 = min(DelV2) ;
Accumulator_size1_actual = max([ ((1+sign(delV_max1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_max1)*perc_threshold).^(1/K)-1) .* delV_max1 , ((1+sign(delV_min1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_min1)*perc_threshold).^(1/K)-1) .* delV_min1 ]) 
Accumulator_size2_actual = max([ ((1+sign(delV_max2)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_max2)*perc_threshold).^(1/K)-1) .* delV_max2 , ((1+sign(delV_min2)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_min2)*perc_threshold).^(1/K)-1) .* delV_min2 ]) 

Total_Accum = Accumulator_size1_actual + Accumulator_size2_actual

figure
plot(t,V,DPt,V_p1,DPt,V_p2), legend('Rail 1 Flow','Rail 2 Flow','Flow delivered to Rail 1 by Pump','Flow delivered to Rail 2 by Pump','Location','NorthWest')
title(['Accumulator Sizes = [',num2str(Accumulator_size1_actual*1000),' , ',num2str(Accumulator_size2_actual*1000) , ']',' Liters'])
ylabel('Volume (m^3)'), xlabel('Time (s)')

%% Plot Fluid volume in each accum
Vf1 = Accumulator_size1_actual - (1-perc_threshold)^(1/K)*Accumulator_size1_actual + DelV1;
Vf2 = Accumulator_size2_actual - (1-perc_threshold)^(1/K)*Accumulator_size2_actual + DelV2;

figure, plot(t,Vf1,t,Vf2)


%% Make table for changing perc_threshold
thresh_vals = [.05 .1 .2] ;
k_vals = [1 1.4];
[a b] = ndgrid(thresh_vals,k_vals);
a = a(:); b = b(:);
for i = 1:length(a)
    perc_threshold = a(i);
    K = b(i);
    Accumulator_Size(i) = max([ ((1+sign(delV_max1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_max1)*perc_threshold).^(1/K)-1) .* delV_max1 , ((1+sign(delV_min1)*perc_threshold)./(1-perc_threshold)).^(1/K) ./ ((1+sign(delV_min1)*perc_threshold).^(1/K)-1) .* delV_min1 ]) ;
end
[a b Accumulator_Size'*1000]



