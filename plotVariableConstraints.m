% This function will plot the inventory and injection/withdrawal curves
% along with the constraints for those values for a piecewise linear model
%
% Inputs:
%   
%   d:  a vector of length n with positive integer values where d(k) for 
%       1 <= k <= n indicated the number of contracts delivered by us in
%       month k (in other words, contracts we sell)
%
%   e:  a vector of length n with positive integer values where d(k) for 
%       1 <= k <= n indicated the number of contracts delivered by us in
%       month k (in other words, contracts we sell)
%
%   piecewiseConstraints: 
%       a cell array containing the injection and withdrawal constraints
%       as three dimensional matrices, indexed by i,j,k where (i,j)
%       representa a coordinate on the monthly piecewise limits of
%       (inventory level, max injection/withdrawal) and k is the index of
%       the current decision variable
%
%   V0: the initial inventory level of the storage
%
%   Vn: the final inventory level of the storage at the end of month n
%
%   L:  a vector of length n indicating the minimal inventory level of gas
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   U:  a vector of length n indicating the maximal inventory level of gas
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   months: a vector where months(k) is the calendar number of the month
%           represented by variable k (for d_k and e_k)
%
%   cap: a scalar representing the maximal inventory capacity in mmbtu
function plotVariableConstraints(d,e,piecewiseConstraints,V0,Vn,L,U,months,cap)

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

i=piecewiseConstraints{1};
w=piecewiseConstraints{2};


%% Calculating inventory levels and monthly "derivative" (injection/...)
for k = 1:length(d)
    v(k) = V0+1e4*(sum(e(1:k-1))-sum(d(1:k-1)));
    injectionMax(k) = i(2,find(i(1,:,k) <= v(k),1, 'Last'),k);
    withdrawalMax(k) = w(2,find(w(1,:,k) >= v(k),1, 'First'),k); 
     
    if(k>1)
        vd(k-1) = v(k)-v(k-1);
    end
end

%% Plotting min/max inventory constraints
figure;
plot(0:length(d), [V0 v], 'LineWidth', 2);

title('Inventory Constraints and Boundary Conditions','fontsize',18)
ylabel('Amount of Natural Gas in Storage (mmbtu)','fontsize',18)
xlabel('Month','fontsize',18);

lowest = min(L(1:length(d)));
highest = max([U(1:length(d)) cap]);

axis([0, length(d), lowest-0.1*highest, 1.1*highest]);
set(gca,'fontsize',14)
hold on;

% Min and max step functions
stairs([0 1:length(d)], [L(1:length(d)-1) Vn Vn], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.5);
plot(1:length(d)-1, L(1:length(d)-1), '.', 'MarkerSize', 20);

stairs([0 1:length(d)], [U(1:length(d)-1) Vn Vn], 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
plot(1:length(d)-1, U(1:length(d)-1), '.', 'MarkerSize', 20);

% Plotting equality constraints
plot(0,V0,'.', 'MarkerSize', 25);
plot(length(d),Vn,'.', 'MarkerSize', 25);
legend({'Inventory','End of Month Inventory Minimum','Inventory Minimum',...
        'End of Month Inventory Maximum','Inventory Maximum','V_0','V_n'},...
        'Location', 'west');

%% Plotting monthly injection/withdrawal constraints
figure;
title('Injection and Withdrawal Constraints', 'fontsize', 18);
set(gca,'fontsize',14)
axis([1, length(d), -1.25*max(withdrawalMax), 1.25*max(injectionMax)]);
xlabel('Month', 'fontsize', 18)
ylabel('Change in Inventory (mmbtu/month)','fontsize', 18);
hold on;

stairs(1:length(d), [vd vd(end)],  'LineWidth', 3);
stairs(1:length(d), [injectionMax], 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
stairs(1:length(d), [-withdrawalMax], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.5);
legend({'dV/dt','Injection maximum','Withdrawal maximum'}, 'Location', 'northeast');

clear;
end