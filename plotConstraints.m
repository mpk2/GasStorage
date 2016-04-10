
% This function will plot the inventory and injection/withdrawal curves
% along with the constraints for those values for a simplified linear model
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
%   I:  a vector of length n representing the daily maximum injection rate
%       in mmbtu for month k (1 <= k <= n) [constant]
%
%   W:  a vector of length n representing the daily maximum withdrawal rate
%       in mmbtu for month k (1 <= k <= n) [constant]
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
function plotConstraints(d,e,i,w,V0,Vn,L,U, months, cap)

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

for k = 1:length(d)
    for j = 1:dpm(k)
        v(sum(dpm(1:k-1))+j) = V0+1e4*(sum(e(1:k-1))-sum(d(1:k-1)));
    end
end

vd(1) = (v(dpm(1))-V0)/dpm(1);

for idx = 1:length(d)
    thisMonth = v(sum(dpm(1:idx)));
    
    if(idx>1)
        monthBefore = v(sum(dpm(months(1:idx-1))));
        vd(idx) = (thisMonth-monthBefore)/dpm(months(idx));
    end
   
end

%% Min and max inventory constraints
figure;
plot(0:sum(dpm(1:length(d))), [V0 v], 'LineWidth', 2);
title('Inventory Constraints and Boundary Conditions','fontsize',18)
ylabel('Amount of Natural Gas in Storage','fontsize',18)
xlabel('Day','fontsize',18);

lowest = min(L(1:length(d)));
highest = max([U(1:length(d)) cap]);

axis([0, length(d), lowest-0.1*highest, 1.1*highest]);

set(gca,'fontsize',14)
hold on;

stairs([0 cumsum(dpm(1:length(d)))], [L(1:length(d)-1) Vn Vn], 'Color', 'red', 'LineStyle', '--');
plot(cumsum(dpm(1:length(d)-1)), L(1:length(d)-1), '.', 'MarkerSize', 20);

stairs([0 1:length(d)], [U(1:length(d)-1) Vn Vn], 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
plot(1:length(d)-1, U(1:length(d)-1), '.', 'MarkerSize', 20);

plot(0,V0,'.', 'MarkerSize', 20);
plot(sum(dpm(1:length(d))),Vn,'.', 'MarkerSize', 20);
legend({'Inventory','End of Month Inventory Minimum','Inventory Minimum','V_0','V_n'}, 'Location', 'southeast');

%% Injection and Withdrawal constraints
figure;
title('Injection and Withdrawal Constraints', 'fontsize', 18);
set(gca,'fontsize',14)


axis([1, sum(dpm(1:length(d))), -1.75*max(w(1:length(d))), 1.75*max(i(1:length(d)))]);
xlabel('Day', 'fontsize', 18)
ylabel('Change in Inventory (mmbtu/day)','fontsize', 18);
hold on;

stairs([1 cumsum(dpm(1:length(d)))], [vd vd(end)], 'LineWidth', 2);
stairs([1 cumsum(dpm(1:length(d)))], [i(1:length(d)) i(length(d))], 'Color', [0 0.5 0], 'LineStyle', '--');
stairs([1 cumsum(dpm(1:length(d)))], [-w(1:length(d)) -w(length(d))], 'Color', 'red', 'LineStyle', '--');
legend({'dV/dt','Injection maximum','Withdrawal maximum'}, 'Location', 'southeast');

clear;
end