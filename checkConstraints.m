% valid (boolean)
%   true - All constraints satisfied
%   false - One or more constraints violated
%
% invalidConstraint (scalar)
%   If valid == false, invalidConstraint is the index in 
%     piecewiseConstraints of the first violated constraint. Otherwise it's
%     an empty array.
%
% y (scalar)
%   If valid == false, y is the function value of the violated constraint,
%   evaluated at the input argument point, x.
%
% x (cell aray of 2x1 matrices)
%   Points in a 2-dimensional plane
%
% piecewiseConstraints (array of constraints for each month)
%   Example: piecewiseConstraints = {[0 1 1 2 2 3;
%                                     0 0 1 1 2 2]};

function [valid, invalidConstraint, y] = checkConstraints(x, piecewiseConstraints)

valid = true;
invalidConstraint = -1;
y = NaN;

for monthIndex = 1:length(piecewiseConstraints)
    % last constraint x-value that x is greater than
    xL = find(  x{monthIndex}(1) >= ...
                piecewiseConstraints{monthIndex}(1,:),1,'Last');
    
    % first constraint x-value that is greater than x
    xR = xL + 1;
    
    % (x,y) points on left and right input point x
    x1 = piecewiseConstraints{monthIndex}(1,xL);
    x2 = piecewiseConstraints{monthIndex}(1,xR);
    y1 = piecewiseConstraints{monthIndex}(2,xL);
    y2 = piecewiseConstraints{monthIndex}(2,xR);
    
    % slope of the line of the constraint directly above or below x
    m = (y2 - y1)/(x2 - x1);
    
    % y-intercept of the line of the constraint directly above or below x
    % y = mx + b <=> b = y - mx
    b = y1 - m*x1; 
    
    % valid if the y-value of the input point, x, is below the constraint
    valid = x{monthIndex}(2) <= m*x{monthIndex}(1) + b;

    if ~valid 
        % return the month for which a constraint is violated
        invalidConstraint = monthIndex; 
        
        % return the function value of the violated constraint
        y = m*x{monthIndex}(1) + b;
        return
    end
end
return
