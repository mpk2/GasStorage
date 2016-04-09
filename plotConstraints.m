function plotConstraints(d,e,i,w,V0,Vn,l)

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

for k = 1:length(d)
    for j = 1:dpm(k)
        v(sum(dpm(1:k-1))+j) = V0+1e4*(sum(e(1:k-1))-sum(d(1:k-1))+j/dpm(k)*(e(k)-d(k)));
    end
end

vd(1) = (v(dpm(1))-V0)/dpm(1);

for idx = 1:length(d)
    thisMonth = v(sum(dpm(1:idx)));
    
    if(idx>1)
        monthBefore = v(sum(dpm(1:idx-1)));
    else 
        monthBefore = V0;
    end
    
    vd(idx) = (thisMonth-monthBefore)/dpm(idx);
end

figure;
plot(0:sum(dpm(1:length(d))), [V0 v], 'LineWidth', 2);
title('Minimum Inventory Constraints and Boundary Conditions','fontsize',18)
ylabel('Amount of Natural Gas in Storage','fontsize',18)
xlabel('Day','fontsize',18);
axis([0, sum(dpm(1:length(d))), -1*(max(v)-min(l(1:length(d)))), max(v)*1.5]);
set(gca,'fontsize',14)
hold on;
stairs([0 cumsum(dpm(1:length(d)))], [l(1:length(d)-1) Vn Vn], 'Color', 'red', 'LineStyle', '--');
plot(cumsum(dpm(1:length(d)-1)), l(1:length(d)-1), '.', 'MarkerSize', 20);
plot(0,V0,'.', 'MarkerSize', 20);
plot(sum(dpm(1:length(d))),Vn,'.', 'MarkerSize', 20);
legend({'Inventory','End of Month Inventory Minimum','Inventory Minimum','V_0','V_n'}, 'Location', 'southeast');

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