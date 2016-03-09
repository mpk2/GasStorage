function monthlyConstraints = dailyToMonthly(dailyPiecewise, cap) 

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];


% Injection

I = dailyPiecewise{1};

dailyInjectionModel = [];
    
for j = 1:length(I(1,:))-1

    % Convert I(1,:) into
    % [ 0   0.2     0.5     0.75;       Starting inventory
    %   7   10.5    4       8   ;       Max days at injection rate
    %   20  15      10      3]          Injection rate
    
    inventoryAdded = (I(1,j+1)-I(1,j))*cap;
    
    % Calculate the number of days it would take at the max injection 
    maxDaysForInventoryAdded = inventoryAdded / I(2,j); 
    dailyInjectionModel(:,j) = [I(1,j); maxDaysForInventoryAdded; I(2,j)];
    
    summedDays = cumsum(dailyInjectionModel(2,:));
    
    if(j>1)
        % Now subtract the irrelevant days based on the starting inventory
        summedDays = summedDays - summedDays(j);
        summedDays(summedDays<0) = 0;
    end
    
    totalDays(i,:) = summedDays;
end

    
for i=1:length(dpm)

    monthlyConstraint = zeros(2,length(I(1,:)));
      
    % Create a monthly piecewise for each percentage
    for j=1:length(I(1,:))-1
        
        startingInventory = I(1,j);
        
        maxInjection = 0;
        
        completeIntervalsIdx = find( totalDays(j,:) <= dpm(i) , 1, 'Last');
        
        % Calculate the total injection from the complete intervals and
        % then add the partial interval amount
        maxIntervalInjection = dailyInjectionModel(2,:) .* dailyInjectionModel(3,:);
        
        partialIntervalInjection = (totalDays(completeIntervalsIdx+1,:)-dpm(i)) ...
                                    .* dailyInjectionModel(3,completeIntervalIdx+1);
                                
        maxInjection = sum(maxIntervalInjection(1:completeIntervalsIdx))... 
                        + partialIntervalInjection;
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,i) = [startingInventory ; maxInjection];
    end
    
end

I = monthlyConstraint;


% Withdrawal

W = dailyPiecewise{2};

dailyWithdrawalModel = [];
    
for j = 1:length(W(1,:))-1

    % Convert W(1,:) into
    %
    % [ 0.2     0.5     0.75    1;      Starting inventory
    %   10.5    4       8       5;      Max days at withdrawal rate
    %   10      11      15      20 ]    Withdrawal rate
    
    inventoryWithdrawn = (W(1,j)-W(1,j+1))*cap;
    
    % Calculate the number of days it would take at the max injection 
    maxDaysForInventoryWithdrawn = inventoryWithdrawn / W(2,j); 
    dailyWithdrawalModel(:,j) = [W(1,j); maxDaysForInventoryWithdrawn; W(2,j)];
    
    summedDays = cumsum(dailyWithdrawalModel(2,:));
    
    if(j>1) 
        % Now subtract the irrelevant days based on the starting inventory
        summedDays = summedDays - summedDays(j);
        summedDays(summedDays<0) = 0;
    end
    
    totalDays(i,:) = summedDays;

end

    
for i=1:length(dpm)

    monthlyConstraint = zeros(2,length(I(1,:)));
      
    % Create a monthly piecewise for each percentage
    for j=1:length(W(1,:))-1
        
        startingInventory = W(1,j);
        maxWithdrawal = 0;
        
        completeIntervalsIdx = find( totalDays(j,:) <= dpm(i) , 1, 'Last');
        
        % Calculate the total withdrawal from the complete intervals and
        % then add the partial interval amount
        maxIntervalWithdrawal = dailyWithdrawalModel(2,:) .* dailyWithdrawalModel(3,:);
        
        partialIntervalWithdrawal = (totalDays(completeIntervalsIdx+1,:)-dpm(i)) ...
                                    .* dailyInjectionModel(3,completeIntervalIdx+1);
                                
        maxWithdrawal = sum(maxIntervalWithdrawal(1:completeIntervalsIdx))... 
                        + partialIntervalWithdrawal;
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,i) = [startingInventory ; maxWithdrawal];
    end
    
end

W = monthlyConstraint;

monthlyConstraints = {I, W};

end