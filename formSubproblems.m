function subproblems = formSubproblems(invalidConstraint, curProb, piecewiseConstraints)

    i = invalidConstraint;
    
    % Find the closest point that does satisfy
    
    % Now cut based on this closest point. Divide the invalid constraint
    % into two parts. One part goes to one subproblem, and the other goes
    % to the second subproblem.
    
    
    
    
    % Relax these problems
    prob1 = reformPiecewise(curProb, piecewiseConstraints1);
    prob2 = reformPiecewise(curProb, piecewiseConstraints2);
    
    
    subproblems = {prob1, prob2};
   


end