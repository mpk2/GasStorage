% This function will determine the optimal number of forward contracts to 
% buy / sell over a given period to optimise profits given constraints on 
% gas storage
%
% Inputs:
%   
%   start: the calendar number of the month to start the period
%
%   finish: the calender number of the month to end the period
%
%   F:  a vector of length n representing the current forward curve for
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
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   U:  a vector of length n indicating the maximal inventory level of gas
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   cap: a scalar representing the maximal inventory capacity in mmbtu
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
function [d,e,fval] = optimizeContracts(start, finish, F, I, W, q, p, c, V0, Vn, L, U, cap)

format short

% Define g to be a constant number of mmbtu per day associated with each
% forward contract
g = 1e4;

% Total interval length is just some modular stuff
n = mod(finish-start+1,12);
n(n==0)=12;

% Want to create an array that cycles through 12
% Do this with mod, going to use 'months' as an index set
months = mod(start:start+n-1,12);
months(months==0)=12;

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
% v(s) is volume at beginning of day s
% length(v) = sum(dpm)

% In terms of decision variables, then (j>=1)
% v(s) = g*( sum( (e(1:k-1)-d(1:k-1)) ) + (e(k)-d(k))*g/dpm(k) - (e(1)-d(1))*g )

A = [];
b = [];

% Inventory minimum:
% -L(k) >= -v(s) where s is defined as above, j=1, k=2:n
for k=2:n
    
    % Preallocate an empty row for this constraint
    A(end+1,:) = zeros(1,2*n);
    
    % We want to limit the sth volume, which is equivalent to limiting
    % contracts
    A(end, 1:k-1) = -g;
    A(end, n+1:n+k-1) = g;
    
    % Limit this to inventory minimum
    b(end+1) = L(k)-V0;
end

% Negate everything since we are inverting the constraint
A = -A;
b = -b;

% Maximum inventory constraints
for k=2:n
   % Preallocate an empty row for this constraint
    A(end+1,:) = zeros(1,2*n);
    
    % We want to limit the sth volume, which is equivalent to limiting
    % contracts
    A(end, 1:k-1) = -g;
    A(end, n+1:n+k-1) = g;
    
    % Limit this to inventory minimum
    b(end+1) = U(k-1)-V0; 
end


% Rate limits:
% Max injection:  v(s+1)-v(s) <= i(v(s)) or g*( e(k)-d(k) ) <= I(v(s))
% Max withdrawal: v(s)-v(s+1) <= w(v(s)) or g*( d(k)-e(k) ) <= W(v(s))
% 
% v(s) = g*( sum( (e(1:k-1)-d(1:k-1)) * dpm(1:k-1) ) + (e(k)-d(k))*j )
% v(s+1) = v(s) + g*( e(k)-d(k) )

% Go through all the months and set injection/withdrawal constraints
for k = 1:n
  
    % Preallocate an empty row for this injection constraint
    A(end+1,:) = zeros(1,2*n);
    
    % Subtract the delivered contracts from the exercised contracts for
    % this month
    A(end,k) = -1;
    A(end,n+k) = 1;
    
    % Constrain to the constant daily injection limit for that month
    b(end+1) = I(k)*dpm(months(k))/g;
    
    
    % Preallocate an empty row for this withdrawal constraint
    A(end+1,:) = zeros(1,2*n);
    
    % Subtract the exercised contracts from the delivered contracts for
    % this month
    A(end,k) = 1;
    A(end,n+k) = -1;
    
    % Constrain to the constant daily withdrawal limit for that month
    b(end+1) = W(k)*dpm(months(k))/g;
    
    
    % Make sure not to withdraw/inject more than possible
    %
    % Preallocate an empty row for this injection constraint
    A(end+1,:) = zeros(1,2*n);
    
    % Subtract the delivered contracts from the exercised contracts for
    % this month
    A(end,k) = -1;
    A(end,n+k) = 1;
    
    % Constrain to the amount in inventory
    b(end+1) = I(k)*dpm(months(k))/g;
    
    
    % Preallocate an empty row for this withdrawal constraint
    A(end+1,:) = zeros(1,2*n);
    
    % Subtract the exercised contracts from the delivered contracts for
    % this month
    A(end,k) = 1;
    A(end,n+k) = -1;
    
    % Constrain to the constant daily withdrawal limit for that month
    b(end+1) = W(k)*dpm(months(k))/g;
    
end


% Set up equality constraints
Aeq = [];
beq = [];

% Boundary conditions:
% Start: v(1) = V0
% End: v(sum(dpm(1:n))) = Vn

% Preallocate an empty row for this constraint
Aeq(end+1,:) = zeros(1,2*n);

% Net out the contracts
Aeq(end, 1:n) = -g;
Aeq(end, n+1:2*n) = g;

% By doing it this way, we ensure that both boundary conditions are
% accounted for
Aeq(end,:);
beq(end+1) = (Vn-V0);

% Optimise, setting the lower bound to all zeros and upper bound to inf
%[x fval] = intlinprog(-c, 1:2*n, A, b, Aeq, beq, [d e], inf*ones(1,2*n));
[x, fval] = linprog(-c, A, b, Aeq, beq, [d e], inf*ones(1,2*n));

% Break up into d and e
x(x<eps) = 0;
d(:) = x(1:end/2);
e(:) = x(end/2+1:end);

fval = -fval;

plotConstraints(d,e,I,W,V0,Vn,L,U, months, cap);

end