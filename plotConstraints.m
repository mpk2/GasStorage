function plotConstraints(d,e,i,w,V0,Vn,l)

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

for k = 1:length(d)
    for j = 1:dpm(k)
        v(sum(dpm(1:k-1))+j) = V0+1e4*(dpm(1:k-1)*(e(1:k-1)-d(1:k-1))+j*(e(k)-d(k)));
    end
end

vd(1) = (v(dpm(1))-V0)/dpm(1);

for idx = 2:length(d)
    vd(idx) = (v(sum(dpm(1:idx)))-v(sum(dpm(1:idx-1))))/dpm(idx);
end

figure;
plot(0:sum(dpm(1:length(d))), [V0 v]);
title('Minimum Inventory Constraints and Boundary Conditions')
ylabel('Amount of Natural Gas in Storage')
xlabel('Day');
axis([0, sum(dpm(1:length(d))), -1*(max(v)-min(l(1:length(d)))), max(v)*1.5]);
hold on;
stairs([0 cumsum(dpm(1:length(d)))], [l(1:length(d)-1) Vn Vn], 'Color', 'red', 'LineStyle', '--');
plot(cumsum(dpm(1:length(d)-1)), l(1:length(d)-1), '.', 'MarkerSize', 20);
plot(0,V0,'.', 'MarkerSize', 20);
plot(sum(dpm(1:length(d))),Vn,'.', 'MarkerSize', 20);
legend({'Inventory Volume','End of Month Inventory Minimum','Inventory Minimum','V_0','V_n'}, 'Location', 'southeast');

figure;
title('Injection and Withdrawal Constraints');
axis([1, sum(dpm(1:length(d))), -1.75*max(w(1:length(d))), 1.75*max(i(1:length(d)))]);
xlabel('Day')
ylabel('Change in Volume (mmbtu/day)');
hold on;

stairs([1 cumsum(dpm(1:length(d)))], [vd vd(end)]);
stairs([1 cumsum(dpm(1:length(d)))], [i(1:length(d)) i(length(d))], 'Color', [0 0.5 0], 'LineStyle', '--');
stairs([1 cumsum(dpm(1:length(d)))], [-w(1:length(d)) -w(length(d))], 'Color', 'red', 'LineStyle', '--');
legend({'dV/dt','Injection maximum','Withdrawal maximum'}, 'Location', 'southeast');

clear;
end