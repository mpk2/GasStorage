function monthlyConstraints = dailyToMonthly(dailyPiecewise, cap) 

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

I = dailyPiecewise{1};
W = dailyPiecewise{2};

% Injection
for i=1:length(I)
    
    % Create a monthly piecewise for each percentage
    for startingInventory = 0:0.01:1
    
        maxInjection = 0;
        monthlyConstraint = zeros(2,dpm(i));
        curInventory = startingInventory;
        
        % Go through the number of days in this month
        for j=1:dpm(i)
            
            % last constraint index that this curInventory is greater than
            % i.e. this is the inventory level %
            initialXindex = find(curInventory > I(1,:)*cap,1,'Last');

            % Add the amount of injection to the currentInventory level
            % and the maximum injection for the month so far
            maxInjection = maxInjection + I(2,initialXindex);
            curInventory = curInventory + I(2,initialXindex);
        end
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,j) = [startingInventory ; maxInjection];
    end
    
end

I = monthlyConstraint;


% Withdrawal
for i=1:length(W)
    
    % Create a monthly piecewise for each percentage
    for startingInventory = 0:0.01:1
    
        maxInjection = 0;
        monthlyConstraint = zeros(2,dpm(i));
        curInventory = startingInventory;
        
        % Go through the number of days in this month
        for j=1:dpm(i)
            
            % last constraint index that this curInventory is greater than
            initialXindex = find(curInventory > I(1,:),1,'Last');

            maxInjection = maxInjection + I(initialXindex);
            curInventory = curInventory + I(initialXindex);
        end
        
        % Add the daily maximum at the current inventory level
        monthlyConstraint(:,j) = [startingInventory ; maxInjection];
    end
    
end

W = monthlyConstraint;

monthlyConstraints = {I, W};

end