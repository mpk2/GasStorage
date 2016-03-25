
% cubProb: struct
% piecewiseConstraints: cell array of 2 x n matrices, where each matrix
%   represents a piecewise Constraint, corresponding to the injection and
%   withdrawal rate
% Output: 2 by 1 cell array of struct, where each struct represents
% a subproblem. 
function subproblems = formSubproblems(curProb, splitPoint, month)

    dpm = [31 28 31 30 31 30 31 31 30 31 30 31];
    
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
    A(end, 1:month-1) = -dpm(1:month-1);
    A(end, n+1:n+month-1) = dpm(1:month-1);
    
    % Add in the current day worth of injection/withdrawal (first day of
    % the month)
    A(end, month) = -1;
    A(end, n+month) = 1;

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