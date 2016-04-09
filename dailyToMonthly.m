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

k=0;

for j = 1:length(injectionInventoryLevel)-1

    % Convert I(1,:) into
    % [ 0   0.2     0.5     0.75;       Starting inventory
    %   7   10.5    4       8   ;       Max days at injection rate
    %   20  15      10      3]          Injection rate
    
    inventoryAdded = (injectionInventoryLevel(j+1)-injectionInventoryLevel(j))*cap;
    
    % Calculate the number of days it would take at the max injection 
    maxDaysForInventoryAdded = inventoryAdded / injectionMaxMMBTU(j); 
    
    dailyInjectionModel(:,j+k) = [injectionInventoryLevel(j); 
                                maxDaysForInventoryAdded-(maxDaysForInventoryAdded>28)*(maxDaysForInventoryAdded-28); 
                                injectionMaxMMBTU(j)];
                            

    subIntervals = 0;
    
    while(maxDaysForInventoryAdded > 28)
        subIntervals = subIntervals + 1;
        maxDaysForInventoryAdded = maxDaysForInventoryAdded - 28;
        
        if(maxDaysForInventoryAdded <= 28)
           
            dailyInjectionModel(:,j+k+1:j+k+subIntervals) = ...
                [   (injectionInventoryLevel(j)+(28*(1:subIntervals)*injectionMaxMMBTU(j))/cap); 
                    28*ones(1,subIntervals-1) maxDaysForInventoryAdded; 
                    injectionMaxMMBTU(j)*ones(1,subIntervals)];
                
            k=k+subIntervals;
        end
    end
                         
end

injectionRate = [dailyInjectionModel(3,:) 0];
injectionInterval = [dailyInjectionModel(2,:) 0];
injectionInventoryLevel = [dailyInjectionModel(1,:) 1];
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
    
I = zeros([2,length(injectionInventoryLevel)+1,n]);

% Go through all of the months
for i=1:n
    
    % Create a monthly piecewise for each percentage discontinuity at the
    % daily level
    for j=1:length(injectionInventoryLevel)
        
        startingInventory = injectionInventoryLevel(j)*cap;
        
        maxInjection = 0;
        
        % Doing in order that months are executed (not calendar order)
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
   
    I(:,:,i) = [monthlyConstraint [cap; 0]];
end


% Withdrawal

W = dailyPiecewise{2};

dailyWithdrawalModel = [];
withdrawalInventoryLevel = W(1,:);
withdrawalMaxMMBTU = W(2,:);    

k=0;

for j = 2:length(withdrawalInventoryLevel)

    % Convert I(1,:) into
    % [ 0   0.2     0.5     0.75;       Starting inventory
    %   7   10.5    4       8   ;       Max days at withdrawal rate
    %   20  15      10      3]          Withdrawal rate
    
    inventoryWithdrawn = (withdrawalInventoryLevel(j)-withdrawalInventoryLevel(j-1))*cap;
    
    % Calculate the number of days it would take at the max withdrawal 
    maxDaysForInventoryWithdrawn = inventoryWithdrawn / withdrawalMaxMMBTU(j-1);                         

    subIntervals = 0;
    
    while(maxDaysForInventoryWithdrawn > 28)
        subIntervals = subIntervals + 1;
        maxDaysForInventoryWithdrawn = maxDaysForInventoryWithdrawn - 28;
        
        if(maxDaysForInventoryWithdrawn <= 28)
           
            dailyWithdrawalModel(:,j+k-1:j+k+subIntervals-2) = ...
                [   withdrawalInventoryLevel(j)-((28*(subIntervals:-1:1)*withdrawalMaxMMBTU(j-1))/cap); 
                    maxDaysForInventoryWithdrawn 28*ones(1,subIntervals-1); 
                    withdrawalMaxMMBTU(j-1)*ones(1,subIntervals)];
                
            k=k+subIntervals;
        end
    end
    
    dailyWithdrawalModel(:,j+k-1) = [withdrawalInventoryLevel(j); 
                                maxDaysForInventoryWithdrawn*(subIntervals==0)+28*(subIntervals>0); 
                                withdrawalMaxMMBTU(j-1)];
    
                            
end

withdrawalRate = dailyWithdrawalModel(3,:) ;
withdrawalInterval = dailyWithdrawalModel(2,:);
withdrawalInventoryLevel = dailyWithdrawalModel(1,:);
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
    
W = zeros([2,length(withdrawalInventoryLevel)+1, n]);

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
        
        % This means that a complete index interval did occur, so calculate
        % how many days of total withdrawal that covers
        if(completeIntervalsIdx < length(totalDaysWithdraw(end-j+1,:))+1)
            daysWithdrawnSoFar = totalDaysWithdraw(j,completeIntervalsIdx);
        else
            daysWithdrawnSoFar = 0;
        end
            
        % Basically, if it hasn't covered ALL intervals, there is partial
        if(completeIntervalsIdx > 1)
            partialIntervalWithdrawal = (dpm(months(i))-daysWithdrawnSoFar) ...
                                        * withdrawalRate(completeIntervalsIdx-1);
        else 
            partialIntervalWithdrawal = 0;
        end
       
       
        if(daysWithdrawnSoFar>0)
            maxWithdrawal = sum(maxIntervalWithdrawal(completeIntervalsIdx:j))... 
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

    W(:,:,i) = [[0;0] monthlyConstraint];
end


monthlyConstraints = {I, W};

end