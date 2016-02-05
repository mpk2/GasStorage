% valid (boolean)
%   true - All constraints satisfied
%   false - One or more constraints violated
%
% invalidConstraint (scalar)
%   If valid == false, invalidConstraint is the index in 
%     piecewiseConstraints of the first violated constraint. Otherwise it's
%     an empty array.
%
% x (cell aray of 2x1 matrices)
%   Points in a 2-dimensional plane
%
% piecewiseConstraints (array of constraints for each month)

function [valid, invalidConstraint] = checkConstraints(x, piecewiseConstraints)
valid = true;
invalidConstraint = [];

for j = 1:length(piecewiseConstraints)
    % first constraint x-value that x is greater than
    xL = find(x{j}(1) > piecewiseConstraints{j}(1,:),1,'Last');
    
    % first constraint x-value that is greater than x
    xR = xL + 1;
    
    % (x,y) points on left and right input point x
    x1 = piecewiseConstraints{j}(1,xL);
    x2 = piecewiseConstraints{j}(1,xR);
    y1 = piecewiseConstraints{j}(2,xL);
    y2 = piecewiseConstraints{j}(2,xR);
    
    % slope of the line of the constraint directly above or below x
    m = (y2 - y1)/(x2 - x1);
    
    % y-intercept of the line of the constraint directly above or below x
    % y = mx + b <=> b = y - mx
    b = y1 - m*x1; 
    
    % valid if the y-value of the input point, x, is below the constraint.
    valid = x{j}(2) <= m*x{j}(1) + b;
    
    if ~valid 
        invalidConstraint = j; 
        return 
    end
end
return