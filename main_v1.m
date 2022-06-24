%% This code is the one used for the FPMC paper

% Fill / empty the rails so that the accumulators sizes are minimized.
% This code should be generalized so that any flows can be given


clear, close all
% Load Drive Cycle
% Wheel Loader Drive cycles
load('JCB5T_C0P_3CPR_Flows.mat') % net zero rail is VERY bumpy
%load('CNH_WL_LC_CPP_Flows.mat') % These two CNH flows are very similar
%load('CNH_WL_LC_CPP_Flows.mat')

%JCB big excavator drive cycles
%load('JCB_bigEX_NinetyDeg_C0P_Flows.mat') % slightly bumpy net zero rail
%load('JCB_bigEX_Grading_CP0_Flows.mat') % Smooth net zero rail
%load('JCB_bigEX_trenching_CPP_Flows.mat')

% Posative flow is flow leaving the accumulator
V1 = cumsum(QR_1)*t(2);
V2 = cumsum(QR_2)*t(2);
V3 = cumsum(QR_3)*t(2);
V = [V2 V3];
%V = [V1 V2 V3];
figure, plot(t,V*1e3)
legend('Rail 1','Rail 2'), ylabel('Volume (L)'), xlabel('Time (s)')


DPdt = .1; % How quickly can we switch between rails?

% Set flow rate and cost to switch
Qave = max(abs(V(end,:)/t(end)));
N = 25;
Q_vals = linspace(2,15,N)*Qave;

chi = .1;
K = 1.4; % ratio of specific heats - note that it is capitilized in the functions to not mess with for loop k

%% Loop over various flow rates
include_graphs = 0;
constrain_choices = '';
%constrain_choices = 'PP'; % This forces both rails to be pumped to
for i = 1:length(Q_vals)
    Q = Q_vals(i);
    accum_sizes(i,:) = main(t,V,DPdt,Q,chi,K,include_graphs,constrain_choices) ;
    disp([num2str(i) ' of ' num2str(length(Q_vals)) ' flow rates have been tested']) ;

end


figure, plot(Q_vals*60e3,sum(accum_sizes,2)*1e3)
ylabel('\textbf{Sum of accumulator sizes (L)}','interpreter','latex'), xlabel('\textbf{Flow Rate (LPM)}','interpreter','latex')
legend(['Ts = ' num2str(DPdt*1e3) ' ms']), grid

save(['Accum_size_norm' num2str(DPdt*1e3) '.mat'],'Q_vals','accum_sizes')



%% Analyse the case with a specific flow rate
[~,best_flow_rate_ind] = min( sum(accum_sizes,2) ) ;
best_Q = Q_vals(best_flow_rate_ind);
disp(['Best flow rate is ' num2str(best_Q*60e3,3) ' LPM'])
include_graphs = 1;
accum_sizes = main(t,V,DPdt,best_Q,chi,K,include_graphs,constrain_choices) ;







%% Functions!!!



function [accum_sizes] = main(t,V,DPdt,Q,chi,K,include_graphs,constrain_choices)

dt = t(2);
dtscale = DPdt/dt; % Step through DP at a different time step than the one given by Drive Cycle
DPt = 0:DPdt:t(end);
fine_flow_vec = Q*dt*(0:(dtscale-1))';

[V_options, nn] = make_V_options(V,Q,DPdt);
% nn is the number of discrete volumes which will be used to make up the
% state space. This is determined from the maximum and minimum volumes for each rail

% Initilize matrices
J   = NaN(prod(nn),length(DPt));
ind = NaN(prod(nn),length(DPt));
for i = 1:size(V,2), V_min{i} = J; V_max{i} = J;  end

% make indexers
indexer = make_indexer(nn,size(V,2));

change_in_ind = make_change_in_index_vector(V_options,indexer,size(V,2),nn,Q,DPdt,constrain_choices) ;


% V min and max at the end of the time considered
for i = 1:size(V,2)
    V_min{i}(:,end) = V_options{i}(indexer{i}) - V(end,i);
    V_max{i}(:,end) = V_options{i}(indexer{i}) - V(end,i);
end

% Cost at the end of the time considered
for iii = 1:prod(nn)
    J(iii,end) = 0;
    for i = 1:size(V,2)
        J(iii,end) = J(iii,end) + get_accum_size(V_min{i}(iii,end),V_max{i}(iii,end),chi,K);
    end
end

% There are (size(V,2)+1)^2-size(V,2) options for rail configurations
% Start at t(end)
for k = 1:length(DPt)-1
    [~,t_ind] = min(abs(t-(t(end)-k*DPdt))); % The time indice for the fine time mesh
    for iii = 1:prod(nn) % Each option in the state space
        for ii = 1:((size(V,2)+1)^2-size(V,2)) % Loop though all the possible decisions at this state and time
            [u(ii), V_min_temp(:,ii), V_max_temp(:,ii)] = calculate_cost(J,t_ind,dtscale,V_min,V_max,V_options,indexer,change_in_ind,V,iii,ii,k,nn,fine_flow_vec,chi,K) ;
        end
        [J(iii,end-k),ind(iii,end-k)] = min(u); % Which choice to make if you find yourself at state iii at time end-k
        [V_min, V_max] = update_V_max_min(V_min_temp,V_max_temp,V_min,V_max,ind(iii,end-k),iii,k); % using this decison, what are the new V_min and V_max?
    end
    %display_percent_done(k,length(DPt));
end

% Get optimal solution back out
% We have to start at zero deltaV in each rail
zero_indices = zeros(1,length(indexer{1})) == 0;
for i = 1:size(V,2)
        zero_indices = zero_indices & (V_options{i}(indexer{i}) == 0) ;
end
zero_ind = find(zero_indices); % The index where all rails are at zero deltaV

% Check that cost matches up with what V_min and V_max are
total_accum_size_calc1 = 0;
for i = 1:size(V,2)
    total_accum_size_calc1 = total_accum_size_calc1 + get_accum_size(V_min{i}(zero_ind,1),V_max{i}(zero_ind,1),chi,K);
end
checks = [];
checks = [checks, abs( total_accum_size_calc1 - min(J(zero_ind,1)) ) / total_accum_size_calc1 < 1e-3 ];

% Get vector of relevant indices
decision_ind(1) = ind(zero_ind,1);
next_index = get_next_index(change_in_ind,decision_ind(1),zero_ind,nn) ;
for k = 2:length(DPt)-1
    decision_ind(k) = ind(next_index,k); % the decision now of where we want to be come next time step
    next_index = get_next_index(change_in_ind,decision_ind(k),next_index,nn) ;
end

% What is the volume that the main pump provides over time?
V_MP = zeros(1,size(V,2)) ; 
for k = 1:length(DPt)-1
    [M,P] = valve_orientation(size(V,2),decision_ind(k)) ;
    V_MP_temp = NaN(dtscale,size(V,2));
    for i = 1:size(V,2)
        V_MP_temp(:,i) = V_MP(end,i) + Q*dt*(1:dtscale)*(P==i) - Q*dt*(1:dtscale)*(M==i) ;
    end
    V_MP = [V_MP; V_MP_temp];
end


deltaV = V_MP - V;
for i = 1:size(V,2)
    checks = [checks, ...
        abs( min(deltaV(:,i)) - V_min{i}(zero_ind,1) ) < 1e-9 ...
        & ...
        abs( max(deltaV(:,i)) - V_max{i}(zero_ind,1) ) < 1e-9 ...
        ] ;
end

if sum(checks) == length(checks)
    disp('All checks were successful')
else
    error(['Check ' num2str(find(~checks)) ' failed. See main function for checks.'])
end

for i = 1:size(V,2)
    accum_sizes(i) = get_accum_size(min(deltaV(:,i)),max(deltaV(:,i)),chi,K) ; %m^3
    disp(['Accumulator volume on rail ' num2str(i) ' is: ' num2str(accum_sizes(i)*1e3,3) ' L']);
end
%
if include_graphs == 1
    make_plots(t,V,V_MP,chi,K);
end
%make_table(delV_min,delV_max);

end





function Vt = get_accum_size(Vmin,Vmax,chi,k)
%%
    Vmin = min([0,Vmin]);
    Vmax = max([0,Vmax]);

    Vt_1 = -Vmin*(1-(1-chi)^(1/k))^(-1); % Total volume required to keep Vmin within the pressure tolerances
    
    % Find Vt 2...
    Vt_2 = ( (1+chi)^(1/k)*(Vmin-Vmax) - Vmin ) /(1 - (1+chi)^(1/k) );

    % Total required accumulator volume:
    Vt = max([Vt_1 Vt_2]);
end



function [J, V_min_temp, V_max_temp] = calculate_cost(J,t_ind,dtscale,V_min,V_max,V_options,indexer,change_in_ind,V,iii,ii,k,nn,fine_flow_vec,chi,K)
%% iii is the index we are currently at
% We want to calculate the cost of going to another index indicated by configuration ii
[M,P] = valve_orientation(size(V,2),ii);
next_index = get_next_index(change_in_ind,ii,iii,nn) ;
if isempty(next_index)
    J = NaN ; % If next_index is empty, this means that the rail we are going to is outside of the state space
    V_min_temp = NaN; V_max_temp = NaN;
elseif length(next_index) >1
    error(['next_index is too long. next_index = ' num2str(next_index)])
else
    V_min_temp = NaN(size(V,2),1); V_max_temp = NaN(size(V,2),1);
    for i = 1:size(V,2)
        V_min_temp(i,1) = V_min{i}(next_index,end-k+1); % Intiliaize with the Max and min volume differences from the index we going to on the next time step
        V_max_temp(i,1) = V_max{i}(next_index,end-k+1);
    end
    for i = 1:size(V,2) % integrate over the finer time mesh
        delV = V_options{i}(indexer{i}(iii)) - fine_flow_vec*(i==M) + fine_flow_vec*(i==P) - V((t_ind):(t_ind+dtscale-1),i) ;
        V_min_temp(i,1) = min([V_min_temp(i,1); delV]) ;
        V_max_temp(i,1) = max([V_max_temp(i,1); delV]) ;
    end
    worst_in_section = 0 ;
    for i = 1:size(V,2)
        worst_in_section = worst_in_section + get_accum_size(V_min_temp(i,1), V_max_temp(i,1) ,chi,K) ;
    end
    J = max([J(next_index,end-k+1) , worst_in_section]) ; % do any points between the two DP time steps require larger accumulators?
end
end



function [M,P] = valve_orientation(n,ii)
%% Called by calculate_cost

% There are n+1 rails
% M = 0 means that tank is motoring
% P = 1 means that the first rail is pumping
if ii <= n + 1
    M = 0;
    P = ii-1;
else
    M = floor((ii-2)/n);
    P = mod((ii-2),n);
    if P >= M
        P = P + 1 ;
    end
end
end


function [V_min, V_max] = update_V_max_min(V_min_temp,V_max_temp,V_min,V_max,decision,iii,k)
%% Called in the main program - inside the for loop which goes over each time step
% This code just puts the V_min_temp from the choice we made into the actual V_min thank we will use later
    for i = 1:length(V_min)
        V_min{i}(iii,end-k) = V_min_temp(i,decision);
        V_max{i}(iii,end-k) = V_max_temp(i,decision);
        % decision tells us which for the options we chose this time step (there are (n+1)^2-n options for n rails)
    end
end



function [change_in_ind] = make_change_in_index_vector(V_options,indexer,n,nn,Q,DPdt,constrain_choices)
%% calculate how to go from one position in state space to another given which option you are choosing
ii_vals = 1:(n+1)^2-n ; 
flag = 0;
iii = 1;
change_in_ind = []; % The difference between where you are and where youre going
while flag == 0
    for ii = ii_vals
        [M,P] = valve_orientation(n,ii) ;
        possible_indeces = ones(1,prod(nn));
        for i = 1:n
            temp = V_options{i}(indexer{i}) - (V_options{i}(indexer{i}(iii))- Q*DPdt*(i==M) + Q*DPdt*(i==P)) ;
            possible_indeces = possible_indeces & abs(temp)<1e-9 ;
        end
        next_index = find(possible_indeces) ;
        if ~isempty(next_index)
            change_in_ind(ii) = iii - next_index ;
            % if we are constraining choices, then adjust this vector so that it points to an infeasible location
            if ~isempty(constrain_choices)
                for i = 1:n
                    if (constrain_choices(i) == 'M' && P == i ) || (constrain_choices(i) == 'P' && M == i )
                        change_in_ind(ii) = prod(nn+1);
                    end
                end
            end
            ii_vals(ii==ii_vals) =[]; % now remove ii from the list of indices we need to try
        end
    end
    % if change_in_ind has all the entries necessary, then we can stop
    if isempty(ii_vals)
        flag = 1;
    end


    % Don't let the loop run forver
    if iii >= length(indexer{1})
        flag = -1;
        error('Error occured in the function "make_change_in_index_vector"')
    end
    iii = iii+ 1;
end
end



function next_index = get_next_index(change_in_ind,ii,iii,nn)
%% Called by calculate_cost
% Given initial index iii, what is the net index if we do action ii
next_index = iii - change_in_ind(ii) ;
if next_index > prod(nn) || next_index <= 0
    next_index = [];
end
end


function display_percent_done(k,kf)
%% Called in the main DP for loop looping over time
% Display what percent done the for loop is - only at key percentages, like 10 percent
if kf > 100
    step = 100/(kf-1); % The finest step that is ever made
    events = [1:1:9, 10:10:100];

    event_occur = find( abs( events - k*100/(kf-1) ) < step/2 );
    if ~ isempty(event_occur)
        disp(['Dynamic programing loop is ' num2str(events(event_occur)) '% complete'])
    end
else
    disp(['Dynamic programing loop is ' num2str(round(k*100/(kf-1))) '% complete'])
end
end



function [V_options, nn] = make_V_options(V,Q,DPdt)
%% Called in main program before main DP loop

% nn is the number of discrete volumes which will be used to make up the
% state space. This is determined from the maximum and minimum volumes for each rail
ss_top = max(V,[],1); % Top of state space
ss_bottom = min(V,[],1); % Bottom of state space
nn_pos = ceil( ss_top/Q/DPdt )     + 2 ; % 1x2 vector  - how many steps would it take to get to the top of the state space?
nn_neg = ceil( -ss_bottom/Q/DPdt ) + 2; % The plus two is to give a little room on the top and bottom of the state space

% Make V_options - which is a matrix of possible values of V_p at each time
for i = 1:size(V,2)
    negative_portion = -Q*DPdt*((nn_neg(i)-1):-1:1);
        
    positive_portion = Q*DPdt*(1:nn_pos(i)-1);
    V_options{i} = [negative_portion 0 positive_portion];
    nn(1,i) = size(V_options{i},2);
end
end



function indexer = make_indexer(nn,n)
%% Makes indexer for V_options - this indexer creates a mesh grid so that each index results in a unique point in state space
lhs = [ '[ temp{1} ' ] ;
rhs = [ 'ndgrid( 1:nn(1) ' ] ;
for i = 2:n
    lhs = [ lhs ', temp{' num2str(i) '} ' ] ;
    rhs = [ rhs ', 1:nn(' num2str(i) ') ' ] ;
end
eval( [lhs '] = ' rhs ') ;'] )
for i = 1:n, indexer{i} = temp{i}(:); end
end



function make_table(delV_min,delV_max)
%% Make table for changing perc_threshold
thresh_vals = [.05 .1 .2] ;
k_vals = [1 1.4];
[a, b] = ndgrid(thresh_vals,k_vals);
a = a(:); b = b(:);
for i = 1:length(a)
    chi = a(i);
    K = b(i);
    Accumulator_Size(i) = get_accum_size(delV_min,delV_max,chi,K);
end
[a b Accumulator_Size'*1000]
end




function make_plots(t,V,V_MP,chi,K)
%%
% make legend for volume plot
legend_str1 = ['legend( '];
legend_str2 = [''];
for i = 1:size(V,2)
    if i == 1 
        legend_str1 = [legend_str1 '"Flow required by rail ' num2str(i) '"'] ;
    else
        legend_str1 = [legend_str1 ' , "Flow required by rail ' num2str(i) '"'] ;
    end
    legend_str2 = [legend_str2 ' , "Flow provided to rail ' num2str(i) '"'] ;
end
figure, plot(t,V*1e3,t,V_MP*1e3), ylabel('Volume (L)'), xlabel('Time (s)'), eval([legend_str1 legend_str2 ')' ])


% Plot Fluid volume in each accum
deltaV = V_MP - V ;
for i = 1:size(V,2)
    V_f(:,i) = deltaV(:,i) - min(deltaV(:,i)) ;
end

figure, plot(t,V_f*1e3), xlabel('Time (s)'), ylabel('Volume (L)')

%% plot Pressure in each accum
for i = 1:size(V,2)
    delV_min = min(deltaV(:,i));
    accum_size = get_accum_size(delV_min,max(deltaV(:,i)),chi,K) ;
    percent_P_error(:,i) = ( (accum_size+delV_min) ./ (accum_size - V_f(:,i)) ).^K - 1;
end
figure, plot(t,percent_P_error), xlabel('Time (s)'), ylabel('Percent Error in Pressure'), hold on
plot([0, t(end)],[chi, chi],'k--',[0, t(end)],[-chi, -chi],'--k'), ylim(1.25*[-chi chi]), hold off
end