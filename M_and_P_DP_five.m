%% This code is the one used for the FPMC paper

% Fill / empty the rails so that the accumulators sizes are minimized.

% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.


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
figure, plot(t,V1*1e3,t,V2*1e3)
legend('Rail 1','Rail 2'), ylabel('Volume (L)'), xlabel('Time (s)')

% Set flow rate and cost to switch
Qave = max(abs(V(end,:)/t(end)));
Q = 7.6826e-04;
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
V1_min = J ;
V1_max = J ;
V2_min = J ;
V2_max = J ;

% Make indexers
a = fliplr(cumsum(ones(1,length(V_MP1)))); b = 2*a(1):-1:a(1)+1;
[x,y] = meshgrid(a,b); X = x(:); Y = y(:);

V1_min(:,end) = -(V(end,1)-V_MP(X));
V1_max(:,end) = -(V(end,1)-V_MP(X));
V2_min(:,end) = -(V(end,2)-V_MP(Y));
V2_max(:,end) = -(V(end,2)-V_MP(Y));

for iii = 1:nn^size(V,2)
    J(iii,end) = get_accum_size(V1_min(iii,end),V1_max(iii,end),perc_threshold,K) + get_accum_size(V2_min(iii,end),V2_max(iii,end),perc_threshold,K) ;
end


%
%PercentDone = 0
for k = 1:length(DPt)-1
    [~,t_ind] = min(abs(t-(t(end)-k*DPdt)));
    
    delV1 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    delV2 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    V1_min(1,end-k) = min([V1_min(1,end-k+1), min(delV1)]);
    V1_max(1,end-k) = max([V1_max(1,end-k+1), max(delV1)]);
    V2_min(1,end-k) = min([V2_min(1,end-k+1), min(delV2)]);
    V2_max(1,end-k) = max([V2_max(1,end-k+1), max(delV2)]);

    J(1,end-k) = get_accum_size(V1_min(1,end-k+1),V1_max(1,end-k+1),perc_threshold,K) + get_accum_size(V2_min(1,end-k+1),V2_max(1,end-k+1),perc_threshold,K);

    ind(1,end-k) = 4;
    for j = 2:nn^size(V,2)
        if V_MP(X(j)) == max(V_MP) % The case where rail 1 should not be filled anymore
            [~,ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),inf,inf,J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            
        elseif V_MP(Y(j)) == min (V_MP) % The case where rail 2 should not be emptied anymore
            [~,ind(j,end-k)] = min([inf,inf,J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            %J(j,end-k) = max( [ J(j,end-k) , abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))) ] );            
            
        else
            [~,ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),J(j-1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-1-length(V_MP1),end-k+1)==1|ind(j-1- length(V_MP1),end-k+1)==3|ind(j-1- length(V_MP1),end-k+1)==4),J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            %J(j,end-k) = max([J(j,end-k) , abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))) ] ); 
            %J(j,end-k) = J(j,end-k) + abs(V(t_ind,1)-V_MP(X(j))) + abs(V(t_ind,2)-V_MP(Y(j))); 
            
        end

        delV1 = -(V(t_ind:t_ind+dtscale-1,1)-V_MP(X(j)));
        delV2 = -(V(t_ind:t_ind+dtscale-1,2)-V_MP(Y(j)));

        if ind(j,end-k) == 1
            V1_min(j,end-k) = min([V1_min(j-1,end-k+1), min(delV1)]);
            V1_max(j,end-k) = max([V1_max(j-1,end-k+1), max(delV1)]);
            V2_min(j,end-k) = min([V2_min(j-1,end-k+1), min(delV2)]);
            V2_max(j,end-k) = max([V2_max(j-1,end-k+1), max(delV2)]);
        elseif ind(j,end-k) == 2
            V1_min(j,end-k) = min([V1_min(j-1-length(V_MP1),end-k+1), min(delV1)]);
            V1_max(j,end-k) = max([V1_max(j-1-length(V_MP1),end-k+1), max(delV1)]);
            V2_min(j,end-k) = min([V2_min(j-1-length(V_MP1),end-k+1), min(delV2)]);
            V2_max(j,end-k) = max([V2_max(j-1-length(V_MP1),end-k+1), max(delV2)]);
        elseif ind(j,end-k) == 3
            V1_min(j,end-k) = min([V1_min(j-length(V_MP1),end-k+1), min(delV1)]);
            V1_max(j,end-k) = max([V1_max(j-length(V_MP1),end-k+1), max(delV1)]);
            V2_min(j,end-k) = min([V2_min(j-length(V_MP1),end-k+1), min(delV2)]);
            V2_max(j,end-k) = max([V2_max(j-length(V_MP1),end-k+1), max(delV2)]);
        else
            V1_min(j,end-k) = min([V1_min(j,end-k+1), min(delV1)]);
            V1_max(j,end-k) = max([V1_max(j,end-k+1), max(delV1)]);
            V2_min(j,end-k) = min([V2_min(j,end-k+1), min(delV2)]);
            V2_max(j,end-k) = max([V2_max(j,end-k+1), max(delV2)]);
        end

        J(j,end-k) = get_accum_size(V1_min(j,end-k),V1_max(j,end-k),perc_threshold,K) + get_accum_size(V2_min(j,end-k),V2_max(j,end-k),perc_threshold,K);


    end
    
    
    
    PercentDone = 100*k/length(DPt);
    if rem(PercentDone,10) < 100/length(DPt)
        PercentDone = (round(PercentDone,-1))
    end
end
PercentDone = 100



%% Get optimal solution back out
% We have to start at zero flow
MinCost = min(J(end,1))
get_accum_size(V1_min(end,1),V1_max(end,1),perc_threshold,K) + get_accum_size(V2_min(end,1),V2_max(end,1),perc_threshold,K)

current_ind = NaN(size(J,2));
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
    DelV1(i) =  V_p1(DPt_ind) - V(i,1);
    DelV2(i)=  V_p2(DPt_ind) - V(i,2);
end
%figure, plot(t,DelV1), ylabel('$\Delta V (m^3)$','interpreter','latex'), xlabel('Time (s)','interpreter','latex')
%figure, plot(t,DelV2), ylabel('$\Delta V (m^3)$','interpreter','latex'), xlabel('Time (s)','interpreter','latex')


delV_max1 = max(DelV1) ; [delV_max1, V1_max(end,1)]
delV_min1 = min(DelV1) ; [delV_min1, V1_min(end,1)]
delV_max2 = max(DelV2) ; [delV_max2, V2_max(end,1)]
delV_min2 = min(DelV2) ; [delV_min2, V2_min(end,1)]
Accumulator_size1_actual = get_accum_size(delV_min1,delV_max1,perc_threshold,K);
Accumulator_size2_actual = get_accum_size(delV_min2,delV_max2,perc_threshold,K);

Total_Accum = Accumulator_size1_actual + Accumulator_size2_actual

figure
plot(t,V*1e3,DPt,V_p1*1e3,DPt,V_p2*1e3), legend('Rail 1 Flow','Rail 2 Flow','Flow delivered to Rail 1 by Pump','Flow delivered to Rail 2 by Pump','Location','NorthWest')
%title(['Accumulator Sizes = [',num2str(Accumulator_size1_actual*1000),' , ',num2str(Accumulator_size2_actual*1000) , ']',' Liters'])
ylabel('Volume (L)'), xlabel('Time (s)')

%% Plot Fluid volume in each accum
Vf1 = DelV1 - delV_min1;
Vf2 = DelV2 - delV_min2;

figure, plot(t,Vf1*1e3,t,Vf2*1e3), xlabel('Time (s)'), ylabel('Volume (L)')

%% plot Pressure in each accum
P1_over_Prail = ( (Accumulator_size1_actual+delV_min1) ./ (Accumulator_size1_actual - Vf1) ).^K;
P2_over_Prail = ( (Accumulator_size2_actual+delV_min2) ./ (Accumulator_size2_actual - Vf2) ).^K;

figure, plot(t,P1_over_Prail-1,t,P2_over_Prail-1), xlabel('Time (s)'), ylabel('$\frac{|P(t) - P_{rail}|}{P_{rail}}$','interpreter','latex','FontSize',22.5), hold on
plot([0, t(end)],[perc_threshold, perc_threshold],'k--',[0, t(end)],[-perc_threshold, -perc_threshold],'--k'), ylim(1.25*[-perc_threshold perc_threshold]), hold off

ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','center')


%% Make table for changing perc_threshold
thresh_vals = [.05 .1 .2] ;
k_vals = [1 1.4];
[a b] = ndgrid(thresh_vals,k_vals);
a = a(:); b = b(:);
for i = 1:length(a)
    perc_threshold = a(i);
    K = b(i);
    Accumulator_Size(i) = get_accum_size(delV_min1,delV_max1,perc_threshold,K);
end
[a b Accumulator_Size'*1000]





function Vt = get_accum_size(Vmin,Vmax,chi,k)
    Vmin = min([0,Vmin]);
    Vmax = max([0,Vmax]);

    Vt_1 = -Vmin*(1-(1-chi)^(1/k))^(-1); % Total volume required to keep Vmin within the pressure tolerances
    
    % Find Vt 2...
    Vt_2 = ( (1+chi)^(1/k)*(Vmin-Vmax) - Vmin ) /(1 - (1+chi)^(1/k) );

    % Total required accumulator volume:
    Vt = max([Vt_1 Vt_2]);
end