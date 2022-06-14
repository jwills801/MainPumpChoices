% This code is what was used for the FPMC paper, along with
% M_and_P_DP_five.m
% Fill / empty the rails so that the accumulators sizes are minimized.

clear, close all
% Load Drive Cycle
load('JCB5T_C0P_3CPR_Flows.mat')



% Start at t(end)
dt = t(2)-t(1);
dtscale = 30; % Step through DP at a different time step than the one given by Drive Cycle
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

    for jjj = 1:nn^size(V,2)
        J(jjj,end) = get_accum_size(V1_min(jjj,end),V1_max(jjj,end),perc_threshold,K) + get_accum_size(V2_min(jjj,end),V2_max(jjj,end),perc_threshold,K) ;
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
                [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),inf,inf,J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);

            elseif V_MP(Y(j)) == min (V_MP) % The case where rail 2 should not be emptied anymore
                [J(j,end-k),ind(j,end-k)] = min([inf,inf,J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);

            else
                [J(j,end-k),ind(j,end-k)] = min([J(j-1,end-k+1)+Cost2Switch*(ind(j,end-k+1)==2|ind(j,end-k+1)==3|ind(j,end-k+1)==4),J(j-1-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-1-length(V_MP1),end-k+1)==1|ind(j-1- length(V_MP1),end-k+1)==3|ind(j-1- length(V_MP1),end-k+1)==4),J(j-length(V_MP1),end-k+1)+Cost2Switch*(ind(j-length(V_MP1),end-k+1)==1|ind(j-length(V_MP1),end-k+1)==2|ind(j-length(V_MP1),end-k+1)==4),J(j,end-k+1)+Cost2Switch*(ind(j,end-k+1)==1|ind(j,end-k+1)==2|ind(j,end-k+1)==3)]);

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




    end




    % Get optimal solution back out
    % We have to start at zero flow
    MinCost = min(J(end,1));
    get_accum_size(V1_min(end,1),V1_max(end,1),perc_threshold,K) + get_accum_size(V2_min(end,1),V2_max(end,1),perc_threshold,K);

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
        DelV1(i) =  V_p1(DPt_ind) - V(i,1);
        DelV2(i)=  V_p2(DPt_ind) - V(i,2);
    end

    delV_max1 = max(DelV1) ;
    delV_min1 = min(DelV1) ;
    delV_max2 = max(DelV2) ;
    delV_min2 = min(DelV2) ;



    Accumulator_sizes(iii,:) = [get_accum_size(delV_min1,delV_max1,perc_threshold,K) , get_accum_size(delV_min2,delV_max2,perc_threshold,K)];
    Accumulatorsizenorm(iii) = sum(Accumulator_sizes(iii,:));

    PercentDone = 100*iii/N;
    if rem(PercentDone,10) < 100/N
        PercentDone = (round(PercentDone,-1))
    end
    clear J ind V_MP1 V_MP
end
PercentDone = 100

[~,winning_ind] = min(Accumulatorsizenorm);
BestFlowRate = Qvals(winning_ind)
AccumulatorSizes = Accumulator_sizes(winning_ind,:)

figure, plot(Qvals,Accumulatorsizenorm), ylabel('Sum accumulator sizes ($m^3$)','interpreter','latex'),
xlabel('Flow Rate ($m^3/s$)','interpreter','latex')
grid

return
%%
close all, clear
load Accum_size_norm_100.mat
load Accum_size_norm_200.mat
load Accum_size_norm_300.mat


figure, plot(Qvals,Accumulatorsizenorm_100,Qvals,Accumulatorsizenorm_200,Qvals,Accumulatorsizenorm_300)
ylabel('\textbf{Sum of accumulator sizes ($m^3$)}','interpreter','latex'), xlabel('\textbf{Flow Rate ($m^3/s$)}','interpreter','latex')
legend('Ts = 100 ms','Ts = 200 ms','Ts = 300 ms'), grid





function Vt = get_accum_size(Vmin,Vmax,chi,k)
    Vmin = min([0,Vmin]);
    Vmax = max([0,Vmax]);

    Vt_1 = -Vmin*(1-(1-chi)^(1/k))^(-1); % Total volume required to keep Vmin within the pressure tolerances
    
    % Find Vt 2...
    Vt_2 = ( (1+chi)^(1/k)*(Vmin-Vmax) - Vmin ) /(1 - (1+chi)^(1/k) );

    % Total required accumulator volume:
    Vt = max([Vt_1 Vt_2]);
end
