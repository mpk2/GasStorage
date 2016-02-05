function [valid, invalidConstraint] = checkConstraints(x, piecewiseConst)

    % Rate limits:
    % Max injection:  v(s+1)-v(s) <= i(v(s)) or g*( e(k)-d(k) ) <= I(v(s))
    % Max withdrawal: v(s)-v(s+1) <= w(v(s)) or g*( d(k)-e(k) ) <= W(v(s))
    % 
    % v(s) = g*( sum( (e(1:k-1)-d(1:k-1)) * dpm(1:k-1) ) + (e(k)-d(k))*j )
    % v(s+1) = v(s) + g*( e(k)-d(k) )
    
    % convert this into Ax <= i(Bx) where A is the summation of contracts and
    % i is the piecewise function that spits out the max injection for a
    % given constraint inventory. Bx gives the inventory on day s.
    
    % If it violates a constraint, mark the constraint and return invalid




end