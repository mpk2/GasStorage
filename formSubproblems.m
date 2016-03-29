
% cubProb: struct
% piecewiseConstraints: cell array of 2 x n matrices, where each matrix
%   represents a piecewise Constraint, corresponding to the injection and
%   withdrawal rate
% Output: 2 by 1 cell array of struct, where each struct represents
% a subproblem. 
function subproblems = formSubproblems(start, finish, curProb, splitPoint, monthIndex)

    dpm = [31 28 31 30 31 30 31 31 30 31 30 31];
    g=1e4;
    
    n = mod(finish-start+1,12);
    n(n==0)=12;

    % Want to create an array that cycles through 12
    % Do this with mod, going to use 'months' as an index set
    months = mod(start:start+n-1,12);
    months(months==0)=12;

    
    bound = splitPoint(1);
    
    n = length(curProb.x0)/2;
    
    % Copy problem
    subProbl = curProb;
    subProbu = curProb;
    
    A = subProbl.Aineq;
    
    % Preallocate an empty row for this constraint
    A(end+1,:) = zeros(1,2*n);
    
    % We want to limit the sth volume, which is equivalent to limiting
    % contracts
    A(end, 1:monthIndex-1) = -g;
    A(end, n+(1:monthIndex-1)) = g;
    
    % Add in the current day worth of injection/withdrawal (first day of
    % the month)
    A(end, monthIndex) = -g/dpm(months(monthIndex));
    A(end, n+monthIndex) = g/dpm(months(monthIndex));

    % Adding A to both subproblems
    subProbl.Aineq = A;
    subProbu.Aineq = -A;
    
    % Adding b to both subproblems
    b = subProbl.bineq;
    
    b(end+1) = bound;
    
    subProbl.bineq = b;
    
    b(end) = -1*b(end);
    
    subProbu.bineq = b;
    
    subproblems = {subProbu, subProbl};
   


end