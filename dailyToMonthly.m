function monthlyConstraints = dailyToMonthly(start, finish, dailyPiecewise, cap) 

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

n = 12-start+1+finish;

% Want to create an array that cycles through 12
% Do this with mod, going to use 'months' as an index set
months = mod(start:start+n-1,12);
months(months==0)=12;

% Injection

I = dailyPiecewise{1};

dailyInjectionModel = [];
    
injectionInventoryLevel = I(1,:);
injectionMaxMMBTU = I(2,:);

for j = 1:length(injectionInventoryLevel)-1

    % Convert I(1,:) into
    % [ 0   0.2     0.5     0.75;       Starting inventory
    %   7   10.5    4       8   ;       Max days at injection rate
    %   20  15      10      3]          Injection rate
    
    inventoryAdded = (injectionInventoryLevel(j+1)-injectionInventoryLevel(j))*cap;
    
    % Calculate the number of days it would take at the max injection 
    maxDaysForInventoryAdded = inventoryAdded / injectionMaxMMBTU(j); 
    dailyInjectionModel(:,j) = [injectionInventoryLevel(j); 
                                maxDaysForInventoryAdded; 
                                injectionMaxMMBTU(j)];
end


totalDays = zeros(length(injectionInventoryLevel)-1);

for j = 1:length(injectionInventoryLevel)
   summedDaysInject = cumsum(dailyInjectionModel(2,:));
    
    if(j>1)
        % Now subtract the irrelevant days based on the starting inventory
        % summedDays(j) is the jth point of the piecewise, so subtracting
        % the cumulative sum of days up to that inventory level
        
        % This should yield a vector which is the number of days injected
        % at each inventory level for this starting inventory level
        summedDaysInject = summedDaysInject - summedDaysInject(j-1);
        summedDaysInject(summedDaysInject<0) = 0;
    end
    
    totalDaysInject(j,:) = summedDaysInject;
end

monthlyConstraint = zeros(2,length(injectionInventoryLevel));
injectionRate = dailyInjectionModel(3,:);
injectionInterval = dailyInjectionModel(2,:);
    
I = zeros([2,length(injectionInventoryLevel),n]);
for i=1:n
      
    % Create a monthly piecewise for each percentage
    for j=1:length(injectionInventoryLevel)
        
        startingInventory = injectionInventoryLevel(j)*cap;
        
        maxInjection = 0;
        
        completeIntervalsIdx = find( totalDaysInject(j,:) <= dpm(months(i)) , 1, 'Last');
       
        if(isempty(completeIntervalsIdx)) 
            completeIntervalsIdx = 0;
        end
        
        % Calculate the total injection from the complete intervals and
        % then add the partial interval amount
       
        maxIntervalInjection = injectionInterval .* injectionRate;
        
        if(completeIntervalsIdx < length(injectionRate))
            if(completeIntervalsIdx ~= 0)
                daysInjectedSoFar = totalDaysInject(j,completeIntervalsIdx);
            else
                daysInjectedSoFar = 0;
            end
            
            partialIntervalInjection = (dpm(months(i))-daysInjectedSoFar) ...
                                        * injectionRate(completeIntervalsIdx+1);
        else
            partialIntervalInjection = 0;
        end
        
        maxInjection = sum(maxIntervalInjection(j:completeIntervalsIdx))... 
                        + partialIntervalInjection;
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,j) = [startingInventory ; maxInjection];
    end
   
    I(:,:,i) = monthlyConstraint;
end


% Withdrawal

W = dailyPiecewise{2};

dailyWithdrawalModel = [];
withdrawalInventoryLevel = W(1,:);
withdrawalMaxMMBTU = W(2,:);    

for j = 1:length(withdrawalInventoryLevel)-1

    % Convert W(1,:) into
    %
    % [ 0.2     0.5     0.75    1;      Starting inventory
    %   10.5    4       8       5;      Max days at withdrawal rate
    %   10      11      15      20 ]    Withdrawal rate
    
    inventoryWithdrawn = (withdrawalInventoryLevel(j+1)-withdrawalInventoryLevel(j))*cap;
    
    % Calculate the number of days it would take at the max injection 
    maxDaysForInventoryWithdrawn = inventoryWithdrawn / withdrawalMaxMMBTU(j); 
    dailyWithdrawalModel(:,j) = [   withdrawalInventoryLevel(j); 
                                    maxDaysForInventoryWithdrawn; 
                                    withdrawalMaxMMBTU(j)];

end


withdrawalRate = dailyWithdrawalModel(3,:);
withdrawalInterval = dailyWithdrawalModel(2,:);
totalDaysWithdraw = [];

for j = 1:length(withdrawalInventoryLevel)
    
   % flip this so that the summed days increase towards 0 (the direction in
   % which withdrawal moves the inventory level)
   summedDaysWithdraw = fliplr(cumsum(fliplr(withdrawalInterval)));
   
    if(j>1) 
        % Now subtract the irrelevant days based on the starting inventory
        summedDaysWithdraw = summedDaysWithdraw - summedDaysWithdraw(end-(j-2));
        summedDaysWithdraw(summedDaysWithdraw<0) = 0;
    end
    
    totalDaysWithdraw = [summedDaysWithdraw; totalDaysWithdraw];
end

monthlyConstraint = zeros(2,length(withdrawalInventoryLevel));
    
W = zeros([2,length(withdrawalInventoryLevel), n]);

for i=1:n
      
    % Create a monthly piecewise for each percentage
    for j=1:length(withdrawalInventoryLevel)
        
        startingInventory = withdrawalInventoryLevel(j)*cap;
        
        maxWithdrawal = 0;
        
        % Find how many complete intervals this month will cover for the
        % jth starting inventory level
        completeIntervalsIdx = find( totalDaysWithdraw(j,:) <= dpm(months(i)) , 1, 'First');
       
        % If there is no interval less than the length of the month wide,
        % then just set the index past the length
        if(isempty(completeIntervalsIdx)) 
            completeIntervalsIdx = length(totalDaysWithdraw(j,:))+1;
        end
        
        % Calculate the total withdrawal from the complete intervals and
        % then add the partial interval amount
        maxIntervalWithdrawal = withdrawalInterval .* withdrawalRate;
        
        if(completeIntervalsIdx < length(totalDaysWithdraw(end-j+1,:))+1)
            daysWithdrawnSoFar = totalDaysWithdraw(j,completeIntervalsIdx);
        else
            daysWithdrawnSoFar = 0;
        end
            
        if(completeIntervalsIdx > 1)
            partialIntervalWithdrawal = (dpm(months(i))-daysWithdrawnSoFar) ...
                                        * withdrawalRate(completeIntervalsIdx-1);
        else 
            partialIntervalWithdrawal = 0;
        end
       
        
        if(daysWithdrawnSoFar>0)
            maxWithdrawal = sum(maxIntervalWithdrawal(completeIntervalsIdx:j-1))... 
                           + partialIntervalWithdrawal;
        else 
            maxWithdrawal = partialIntervalWithdrawal;
        end
        
        
        % Just to make sure that you can't withdraw more than you have...
        if(maxWithdrawal > withdrawalInventoryLevel(j)*cap)
            maxWithdrawal = withdrawalInventoryLevel(j)*cap;
        end
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,j) = [startingInventory ; maxWithdrawal];
    end

    W(:,:,i) = monthlyConstraint;
end


monthlyConstraints = {I, W};

end