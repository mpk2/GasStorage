% This function will determine the optimal number of forward contracts to 
% buy / sell over a given period to optimise profits given constraints on 
% gas storage
%
% Inputs:
%   
%   n:  the number of months over which to optimise the contracts
%
%   f:  a vector of length n representing the current forward curve for
%       each month k such that f(k) = forward contract price for month k
%       1 <= k <= n
%
%   I:  a vector of length n representing the daily maximum injection rate
%       in mmbtu for month k (1 <= k <= n) [constant]
%
%   W:  a vector of length n representing the daily maximum withdrawal rate
%       in mmbtu for month k (1 <= k <= n) [constant]
%
%   q:  a function representing the price-dependent cost of withdrawal
%
%   p:  a function representing the price-dependent cost of injection
%
%   c:  a vector of length n representing the month-dependent cost of
%       injection/withdrawal
%       
%   V0: the initial inventory level of the storage
%
%   Vn: the final inventory level of the storage at the end of month n
%
%   L:  a vector of length n indicating the minimal inventory level of gas
%       required to be kept during month k (1 <= k <= n)
%
% Outputs:
%
%   d:  a vector of length n with positive integer values where d(k) for 
%       1 <= k <= n indicated the number of contracts delivered by us in
%       month k (in other words, contracts we sell)
%
%   e:  a vector of length n with positive integer values where d(k) for 
%       1 <= k <= n indicated the number of contracts delivered by us in
%       month k (in other words, contracts we sell)
function [d,e] = optimizeContracts(n, F, I, W, q, p, c, V0, Vn, L)

% Define g to be a constant number of mmbtu per day associated with each
% forward contract
g = 1e5;

% Number of days in each month!
dpm = [31 28 31 30 31 30 31 31 30 31 30 31];

% Decision varibles of length n
% d(k) indicates number of foward contracts delivered during month k
% e(k) indicated number of forward contracts exercised during month k
% d implies future withdrawal from our gas supply
% e implies future replenishment of our gas supply
d = zeros(n,1);
e = zeros(n,1);

% Set P and Q to be cost vectors including c
P = p(F) + c;
Q = q(F) + c;

% Objective vector -- this makes it so we optimise F*(d-e)-g*(Q*d-P*e)
% In other words, we optimise (gas revenue - gas expenditure - gas
% injection/withdrawal costs)
c = [F-g*Q; -F-g*P];

% This is in general helpful for calculations, but it's more helpful
% perhaps when we need to move into piecewise
% 
% s=sum(dpm(1:k-1))+j; k=month, j=day
% s is total number of days that have passed
%
% In order to go from s->k,j
% cumsum(dpm) and find first where >= s, take k=index-1
% then j=s-sum(dpm(1:k))
%
% v(s) is volume for day s
% length(v) = sum(dpm)

% In terms of decision variables, then (j>=1)
% v(s) = v(1) + g*( sum( (e(1:k-1)-d(1:k-1)) * dpm(1:k-1) ) + (e(k)-d(k))*j - (e(1)-d(1)) )
% v(1) = V0

A = [];
b = [];

% Inventory minimum:
% -L(k) >= -v(s) where s is defined as above, j=1, k=1:n
for k=1:n
    
    % Preallocate an empty row for this constraint
    A(end+1,:) = zeros(1,n);
    
    % We want to limit the sth volume, which is equivalent to limiting
    % contracts
    A(end, 1:k-1) = dpm(1:k-1);
    A(end, n+1:2*k-1) = -dpm(1:k-1);
    
    % Add in the current day worth of injection/withdrawal (first day of
    % the month)
    A(end, k) = 1;
    A(end, 2*k) = 1;
    
    % Make sure not to double count first day of year (s=1)
    A(end, 1) = A(end,1) - 1;
    A(end, n+1) = A(end,1) - 1;
    
    % Negate everything since we are inverting the constraint
    A(end,:) = - A(end,:);
    
    % Limit this to inventory minimum
    b(end+1) = -(L(k)-V0)/g;
end



% Rate limits:
% Max injection:  v(s+1)-v(s) <= i(v(s)) or g*( e(k)-d(k) ) <= I(v(s))
% Max withdrawal: v(s)-v(s+1) <= w(v(s)) or g*( d(k)-e(k) ) <= W(v(s))
% 
% v(s) = v(1) + g*( sum( (e(1:k-1)-d(1:k-1)) * dpm(1:k-1) ) + (e(k)-d(k))*j - (e(1)-d(1)) )
% v(s+1) = v(s) + g*( e(k)-d(k) )

% Go through all the months and set injection/withdrawal constraints
for k = 1:n
  
    % Preallocate an empty row for this injection constraint
    A(end+1,:) = zeros(1,n);
    
    % Subtract the delivered contracts from the exercised contracts for
    % this month
    A(end,k) = -1;
    A(end,n+k) = 1;
    
    % Constrain to the constant daily injection limit for that month
    b(end+1) = I(k);
    
    
    % Preallocate an empty row for this withdrawal constraint
    A(end+1,:) = zeros(1,n);
    
    % Subtract the exercised contracts from the delivered contracts for
    % this month
    A(end,k) = 1;
    A(end,n+k) = -1;
    
    % Constrain to the constant daily withdrawal limit for that month
    b(end+1) = W(k);
end

% Set up equality constraints
Aeq = [];
beq = [];

% Boundary conditions:
% Start: v(1) = V0
% End: v(sum(dpm(1:n))) = Vn

% The volume v(1) for k=1, j=1 should be V0
% V0 = v(1) + g*((e(1)-d(1))*1 - (e(1)-d(1)) ) = v(1)
Aeq(end+1,:) = zeros(1:n);
beq(end+1) = V0; 

% Preallocate an empty row for this constraint
Aeq(end+1,:) = zeros(1,n);

% We want to set the last volume to Vn, so add up all contracts in first
% previous full months
Aeq(end, 1:n-1) = dpm(1:k-1);
Aeq(end, n+1:2*n-1) = -dpm(1:k-1);

% Add in the total days in the month worth of injection/withdrawal
Aeq(end, n) = dpm(n);
Aeq(end, 2*n) = dpm(n);

% Make sure not to double count first day of year (s=1)
Aeq(end, 1) = Aeq(end,1) - 1;
Aeq(end, n+1) = Aeq(end,1) - 1;

% Limit this to inventory maximum
beq(end+1) = (Vn-V0)/g;

% Optimise, setting the lower bound to all zeros and upper bound to inf
[x, fval] = linprog(c, A, b, Aeq, beq, [d e], inf*ones(1,n));

% Break up into d and e
d(:) = x(1:n);
e(:) = x(n+1:2*n);

% Display the solution vector and the evaluation for debugging 
disp(x);
disp(fval);

end