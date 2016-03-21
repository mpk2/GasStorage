function convexProblem = reformPiecewise(initProb, piecewiseConstraints)

convexProblem = initProb;
n = length(initProb.lb)/2;

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

for month = 1:n

    for constraint=1
        
        x = [0 piecewiseConstraints{constraint}(1,:,month)];
        y = [0 piecewiseConstraints{constraint}(2,:,month)];

        % piecewise constraints are of form [x;y]
        k = convhull(x, y);
        
        
        % Go through all the segments counterclockwise, except the first one (0)
        for i=length(k)-(2-constraint):-1:2

            % Extract (x1,y1) (x2,y2) from the convex hull info
            x1 = x(k(i));
            y1 = y(k(i));
            x2 = x(k(i-1));
            y2 = y(k(i-1));

            % Form the individual linear constraints
            m = (y1-y2) / (x1-x2); % -1
            b = x1-m*y2; % 0

            % Preallocate an empty row for this constraint
            inventoryVector = zeros(1,2*n);

            % We want to limit the sth volume, which is equivalent to limiting
            % contracts
            inventoryVector(1:month-1) = -dpm(1:month-1);
            inventoryVector(n+1:n+month-1) = dpm(1:month-1);

            % Add in the current day worth of injection/withdrawal (first day of
            % the month)
            inventoryVector(month) = -1;
            inventoryVector(n+month) = 1;

            % Preallocate an empty row for this injection/withdrawal constraint
            deltaVector = zeros(1,2*n);

            % Subtract the delivered contracts from the exercised contracts for
            % this month for injection, reverse for withdrawal
            deltaVector(month) = -1^(constraint-1);
            deltaVector(n+month) = -1^(2*constraint-1);

            newA = (-1)^constraint*(deltaVector+m*inventoryVector);

            convexProblem.Aineq = [convexProblem.Aineq; newA];
            convexProblem.bineq = [convexProblem.bineq, (-1)^constraint*-b];
        end

    end

end
