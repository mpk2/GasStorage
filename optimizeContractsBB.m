% This function will determine the optimal number of forward contracts to 
% buy / sell over a given period to optimise profits given constraints on 
% gas storage
%
% Inputs:
%   
%   start/finish:  calendar 'number' of the start and finish months
%
%   F:  a vector of length n representing the current forward curve for
%       each month k such that f(k) = forward contract price for month k
%       1 <= k <= n
%
%   I:  an array containing pairs of vectors representing ordered pairs 
%       that define boundary points in the daily maximum injection rate 
%       function as f(inventory level) = max injection rate in mmbtu
%
%   W: an array containing pairs of vectors representing ordered pairs
%      that define boundary points in the daily maximum injection rate
%      function as f(inventory level) = max withdrawal rate in mmbtu
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
%   cap: a scalar representing the maximal inventory capacity
%
function [d, e, fval] = optimizeContractsBB(start, finish, F, I, W, q, p, c, V0, Vn, L, cap)

g=1e4;

% Develop the original problem (without injection or withdrawal
% constraints)
initProb = formProblem(start, finish, F, q, p, c, V0, Vn, L);

% Save the piecewise constraints
dailyPiecewiseConstraints = {I, W};
c = initProb.f;

% Turn these into monthly constraints...
piecewiseConstraints = dailyToMonthly(dailyPiecewiseConstraints, cap);

% Form the convex hull of the constraints
relaxedProb = reformPiecewise(initProb, piecewiseConstraints);

% Begin the stack
LIST = [relaxedProb];

% Total interval length
n = 12-start+1+finish;

% Initialise upper bound based on x=zeros
x = zeros(2*n,1);

curOptimal = c*x;

% Pop off the stack until it's empty
while (~isempty(LIST))
    
    % Pull off the first problem 
    curProblem = LIST(:,1);
    LIST(:,1) = [];
    
    % Calculate the optimisation to this problem
    [x_s,~,flag] = linprog(curProblem);
    
    % if it cannot be pruned by infeasibility or bound (i.e. is lower than
    % the current best legitimate candidate)
    if (~isempty(x_s) && c*x_s < curOptimal && flag == 1)
        
        % Need to create cell array of inventory vs injection for each
        % month
        d_s = x_s(1:end/2);
        e_s = x_s(end/2+1:end);
        delta = [(V0+cumsum(e_s-d_s)*g) (e_s-d_s)*g];
        delta(:,1) = delta(:,1) - delta(:,2);
        delta = mat2cell(delta,ones(length(delta),1),2);
        
        relevantConstraints = {};
        
        for month=1:size(delta,1)
           relevantConstraints{month} = piecewiseConstraints{1 + (delta{month}(2) > 0)}(:,:,month);
           delta{month}(2) = abs(delta{month}(2)); 
        end 
        
        % Check against piecewise constraints
        [valid, invalidConstraint] = checkConstraints(delta, relevantConstraints);
        
        % If it satisfied the constraints (and is greater from before)
        if(valid)
            curOptimal = c*x_s;
            x = x_s;
            
        % It didn't satisfy constraints and is still greater, branch
        else
            
            % Subdivide the problem into two on either side of the point based
            % on the constraint that was violated (I or W, then which
            % segment)
            subProblems = formSubproblems(relaxedProb, delta{invalidConstraint}/cap, ...
                                          invalidConstraint);
            
            % breadth first
            LIST = [LIST subProblems{1} subProblems{2}];
            
        end
    end
end

fval = curOptimal;
 d = x(1:end/2);
 e = x(end/2+1:end);
return