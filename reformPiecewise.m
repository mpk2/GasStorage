function convexProblem = reformPiecewise(start, finish, cap, V0, initProb, piecewiseConstraints)

convexProblem = initProb;
n = mod(finish-start+1,12);
n(n==0)=12;

% Want to create an array that cycles through 12
% Do this with mod, going to use 'months' as an index set
months = mod(start:start+n-1,12);
months(months==0)=12;


g=1e4;

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

for monthIndex = 1:n

    for constraint=1:2
        
        x = [0 piecewiseConstraints{constraint}(1,:,monthIndex)];
        y = [-1 piecewiseConstraints{constraint}(2,:,monthIndex)];

        
        % piecewise constraints are of form [x;y]
        k = convhull(x, y);
        k = k(k~=1);
        
%        figure
%         plot(x(k(length(k):-1:1)),...
%         y(k(length(k):-1:1)),'r-',x,y,'b*')
        
        %figure
        % Go through all the segments counterclockwise, except the first one (i=1)
        for i=length(k):-1:2

            % Extract (x1,y1) (x2,y2) from the convex hull info
            x1 = x(k(i-1));
            y1 = y(k(i-1));
            x2 = x(k(i));
            y2 = y(k(i));

            % Form the individual linear constraints
            m = (y2-y1) / (x2-x1);
            b = -x1*m+y1;

            %xx = linspace(0,1e6,1e6);
            %plot(xx,m*xx+b);
            %ylim([0,4e5])
            %hold on;
            % Preallocate an empty row for this constraint
            v = zeros(1,2*n);

            % Basically we have -w(inventory) <= dInventory <= i(inventory)
            % dInventory = delta*X
            % inventory = v*X + V0
            % i(inventory) = m(v*X+V0) + b [same form for w(inventory)]
            v(1:monthIndex-1) = -g;
            v(n+(1:monthIndex-1)) = g;

            % Add in the current day worth of injection/withdrawal (first day of
            % the month)
            v(monthIndex) = -g/dpm(months(monthIndex));
            v(n+monthIndex) = g/dpm(months(monthIndex));

            % Preallocate an empty row for this injection/withdrawal constraint
            delta = zeros(1,2*n);

            % Subtract the delivered contracts from the exercised contracts for
            % this month for injection, reverse for withdrawal
            delta(monthIndex) = -g;
            delta(n+monthIndex) = g;

            newA = (-1)^(constraint-1)*(delta)-m*v;

            convexProblem.Aineq = [convexProblem.Aineq; newA];
            convexProblem.bineq = [convexProblem.bineq, b+(m*V0)];
        end
        

    end
    
end


end
