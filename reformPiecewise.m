function convexProblem = reformPiecewise(initProb, piecewiseConstraints)

x = [0 piecewiseConstraints{1}(1,:)];
y = [0 piecewiseConstraints{1}(2,:)];


convexProblem = initProb;

% piecewise constraints are of form [x;y]
k = convhull(x, y);

n = length(initProb.lb)/2;

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

convexProb = initProb;

% Go through all the segments counterclockwise, except the first one (0)
for i=length(k):-1:3
    
    % Extract (x1,y1) (x2,y2) from the convex hull info
    x1 = x(k(i));
    y1 = y(k(i));
    x2 = x(k(i-1));
    y2 = y(k(i-1));
    
    % Form the individual linear constraints
    m = (y1-y2) / (x1-x2);
    b = x1-m*y2;
    
    % Over each month
    for j = 1:n
        
        % Preallocate an empty row for this constraint
        inventoryVector = zeros(1,2*n);

        % We want to limit the sth volume, which is equivalent to limiting
        % contracts
        inventoryVector(1:j-1) = -dpm(1:j-1);
        inventoryVector(n+1:n+j-1) = dpm(1:j-1);

        % Add in the current day worth of injection/withdrawal (first day of
        % the month)
        inventoryVector(k) = -1;
        inventoryVector(n+j) = 1;

        % Preallocate an empty row for this injection constraint
        injectionVector = zeros(1,2*n);
        withdrawalVector = zeros(1,2*n);

        % Subtract the delivered contracts from the exercised contracts for
        % this month
        injectionVector(j) = -1;
        injectionVector(n+j) = 1;

        % Subtract the exercised contracts from the delivered contracts for
        % this month
        withdrawalVector(j) = 1;
        withdrawalVector(n+j) = -1;
        
        newA = injectionVector-m*inventoryVector;
        newA(end,:) = -withdrawalVector + m*inventoryVector;
        
        convexProb.Aineq = [convexProb.Aineq; newA];
        convexProb.bineq = [convexProb.bineq, b, -b];
    end

end

end