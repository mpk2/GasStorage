%
%   piecewiseProblem has properties
%
%
function convexProblem = reformPiecewise(initProb, piecewiseConstraints)

% piecewise constraints are of form [x;y]
k = convhull([0 piecewiseConstraints(1,:)], [0 piecewiseConstraints(2,:)]);

% Go through all the segments
for i=1:length(k)
    
    % Extract (x1,y1) (x2,y2) from the convex hull info
    
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
    end

    initProb.Aineq = [initProb.Aineq; newA];
    initProb.b = [initProb.b b];
end

end