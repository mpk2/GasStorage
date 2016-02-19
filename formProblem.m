function gasProblem = formProblem(n, F, q, p, c, V0, Vn, L)

% Define g to be a constant number of mmbtu per day associated with each
% forward contract
g = 1e4;

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
c = [F-g*Q -F-g*P];

A = [];
b = [];

% Inventory minimum:
% -L(k) >= -v(s) where s is defined as above, j=1, k=2:n
for k=2:n
    % Preallocate an empty row for this constraint
    A(end+1,:) = zeros(1,2*n);
    
    % We want to limit the sth volume, which is equivalent to limiting
    % contracts
    A(end, 1:k-1) = -dpm(1:k-1);
    A(end, n+1:n+k-1) = dpm(1:k-1);
    
    % Add in the current day worth of injection/withdrawal (first day of
    % the month)
    A(end, k) = -1;
    A(end, n+k) = 1;
    
    % Limit this to inventory minimum
    b(end+1) = (L(k-1)-V0)/g;
end

% Negate everything since we are inverting the constraint
A = -A;
b = -b;

% Set up equality constraints
Aeq = [];
beq = [];

% Boundary conditions:
% Start: v(1) = V0
% End: v(sum(dpm(1:n))) = Vn

% Preallocate an empty row for this constraint
Aeq(end+1,:) = zeros(1,2*n);

% Add up all net contract-days
Aeq(end, 1:n) = -dpm(1:n);
Aeq(end, n+1:2*n) = dpm(1:n);

% By doing it this way, we ensure that both boundary conditions are
% accounted for
beq(end+1) = (Vn-V0)/g;


% Building Prob struct
options = optimoptions('linprog','Display','off');
gasProblem = struct('x0',zeros(1, 2*n),'Aeq',Aeq,'beq',beq,...
    'f',-c','Aineq',A,'bineq',b,'lb',zeros(1, 2*n),'ub',inf*ones(1,2*n),...
    'solver','linprog','options',options);

end