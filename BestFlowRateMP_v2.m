
% Fill / empty the rails so that the accumulators sizes are minimized.

% This version has adjusted the Cost2Switch to hopefully do its job better,
% although the success of this feature has not been tested and is hard to
% see in regular drive cycles. A Drive cycle could be cooked up to test
% whether switches are correctly penalyzed.

% This code takes a drive cycle, and modifies it so that it has two rails
% that need main pump flow

%clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')



% Start at t(end)
dt = t(2)-t(1);
dtscale = 10; % Step through DP at a different time step than the one given by Drive Cycle
DPdt = dtscale*dt;
DPt = 0:DPdt:t(end);

perc_threshold = .1;
K = 1.4;

% figure(1), plot(t,cumsum(QR_1)*dt,t,cumsum(QR_2)*dt,t,cumsum(QR_3)*dt)
% legend('Rail 1','Rail 2','Rail 3'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')

% Posative flow is flow leaving the accumulator

% We are pumping to rail 3 and inverting rail 1 and pumping to it
V1 = cumsum(QR_3)*dt;
V2 = cumsum(QR_1)*dt;
V = [V1 V2];
% figure(1), plot(t,V1,t,V2)
% legend('Rail 1','Rail 2'), ylabel('Cummulative Flow (m^3)'), xlabel('Time (s)')

% Set flow rate and cost to switch
Qave = max(abs(V(end,:)/t(end)));
N = 100;
Qvals = linspace(1*Qave,20*Qave,N);
Cost2Switch = 0*1e-2;

PercentDone = 0
for iii = 1:N
    
    
% nn is the number of time steps with flow it takes to get the the required
% volume
nn = max(ceil(V1(end)/Qvals(iii)/DPdt) +2,ceil(V2(end)/-Qvals(iii)/DPdt) +2);


% Make V_MP - which is a matrix of possible values of V_p at each time
for i = 0:nn-1
    V_MP1(i+1) = i*Qvals(iii)*DPdt;
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
    disp('Size of V not accounted for in this code')
    return
elseif size(V,2) == 2
    [x,y] = meshgrid(a,b); X = x(:); Y = y(:);

    delV1 = -(V(end,1)-V_MP(X));
    delV2 = -(V(end,2)-V_MP(Y));
    Accumulator_size1(:,end) = -2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV1.*(delV1>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV1.*(delV1<0);
    Accumulator_size2(:,end) = -2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV2.*(delV2>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV2.*(delV2<0);

    J(:,end) = Accumulator_size1(:,end) + Accumulator_size2(:,end) ;
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
    
    
    delV1 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    delV2 = -(V1(t_ind:t_ind+dtscale-1)-V_MP(1));
    A1 = max(-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV1.*(delV1>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV1.*(delV1<0));
    A2 = max(-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV2.*(delV2>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV2.*(delV2<0));
    Accumulator_size1(1,end-k) = max([A1 Accumulator_size1(1,end-k+1)]) ;
    Accumulator_size2(1,end-k) = max([A2 Accumulator_size2(1,end-k+1)]) ;
    J(1,end-k) = Accumulator_size1(1,end-k+1) + Accumulator_size2(1,end-k+1);

    ind(1,end-k) = 4;
    for j = 2:nn^size(V,2)
        if V_MP(X(j)) == max(V_MP) % The case where rail 1 should not be filled anymore
            [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),inf,inf,J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            
        elseif V_MP(Y(j)) == min (V_MP) % The case where rail 2 should not be emptied anymore
            [J(j,end-k),ind(j,end-k)] = min([inf,inf,J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            
        else
            [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),J(j-1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-1-length(V_MP1),end-k+1)==1|ind(j-1- length(V_MP1),end-k+1)==3|ind(j-1- length(V_MP1),end-k+1)==4),J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);
            
        end

        delV1 = -(V(t_ind:t_ind+dtscale-1,1)-V_MP(X(j)));
        delV2 = -(V(t_ind:t_ind+dtscale-1,2)-V_MP(Y(j)));
        A1 = max(-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV1.*(delV1>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV1.*(delV1<0));
        A2 = max(-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV2.*(delV2>0) -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV2.*(delV2<0));
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
    
    
    

end




% Get optimal solution back out
% We have to start at zero flow
MinCost= min(J(end,1));
Accumulator_size1(end,1) + Accumulator_size2(end,1) ;

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
    V_p1(jj+1) = V_p1(jj) + (ind(current_ind(jj),jj)==2|ind(current_ind(jj),jj)==3)*Qvals(iii)*DPdt;
    V_p2(jj+1) = V_p2(jj) - (ind(current_ind(jj),jj)==1|ind(current_ind(jj),jj)==2)*Qvals(iii)*DPdt;
end




% What is the required accumulator size?
for i = 1:length(t)
    [~,DPt_ind] = min(abs( DPt - t(i) + DPdt/2 - t(2)));
    Error1(i) = V(i,1) - V_p1(DPt_ind);
    Error2(i) = V(i,2) - V_p2(DPt_ind);
end
DelV1 = -Error1 ;
DelV2 = -Error2 ;

delV_max1 = max(DelV1) ;
delV_min1 = min(DelV1) ;
delV_max2 = max(DelV2) ;
delV_min2 = min(DelV2) ;
Accumulator_size1_actual = max([-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV_max1 , -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV_min1]) ;
Accumulator_size2_actual = max([-2*(1+perc_threshold)^(1/K)/(1-(1+perc_threshold)^(1/K))*delV_max2 , -2*(1-perc_threshold)^(1/K)/(1-(1-perc_threshold)^(1/K))*delV_min2]) ;


Accumulator_sizes(iii,:) = [Accumulator_size1_actual Accumulator_size2_actual];
Accumulatorsizenorm(iii) = norm([Accumulator_size1_actual Accumulator_size2_actual],1);

    PercentDone = 100*iii/N;
    if rem(PercentDone,10) < 100/N
        PercentDone = (round(PercentDone,-1))
    end
    clear J ind V_MP1 V_MP Accumulator_size1 Accumulator_size2
end
PercentDone = 100

[~,winning_ind] = min(Accumulatorsizenorm);
BestFlowRate = Qvals(winning_ind)
AccumulatorSizes = Accumulator_sizes(winning_ind,:)

figure, plot(Qvals,Accumulatorsizenorm), ylabel('1 norm of accumulator sizes ($m^3$)','interpreter','latex'),
xlabel('Flow Rate ($m^3/s$)','interpreter','latex')
grid

return
%%
close all, clear
load Accum_size_norm_100.mat
load Accum_size_norm_200.mat
load Accum_size_norm_300.mat


figure, plot(Qvals,Accumulatorsizenorm_100,Qvals,Accumulatorsizenorm_200,Qvals,Accumulatorsizenorm_300)
ylabel('\textbf{1 norm of accumulator sizes ($m^3$)}','interpreter','latex'), xlabel('\textbf{Flow Rate ($m^3/s$)}','interpreter','latex')
legend('Ts = 100 ms','Ts = 200 ms','Ts = 300 ms'), grid


