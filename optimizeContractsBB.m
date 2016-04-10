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
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   U:  a vector of length n indicating the maximal inventory level of gas
%       required to be present at the end of the last day of  month k (1 <= k <= n)
%
%   cap: a scalar representing the maximal inventory capacity
%
function [d, e, fval] = optimizeContractsBB(start, finish, F, I, W, q, p, c, ...
                                            V0, Vn, L, U, cap)

g=1e4;

% Develop the original problem (without injection or withdrawal
% constraints)
initProb = formProblem(start, finish, F, q, p, c, V0, Vn, L, U);

% Save the piecewise constraints
dailyPiecewiseConstraints = {I, W};
c = initProb.f;

% Turn these into monthly constraints...
piecewiseConstraints = dailyToMonthly(start, finish, dailyPiecewiseConstraints, cap);

% Form the convex hull of the constraints
relaxedProb = reformPiecewise(start, finish, cap, V0, initProb, piecewiseConstraints);

% Begin the stack
STACK = [relaxedProb];

% Total interval length is just some modular stuff
n = mod(finish-start+1,12);
n(n==0)=12;

% Want to create an array that cycles through 12
% Do this with mod, going to use 'months' as an index set
months = mod(start:start+n-1,12);
months(months==0)=12;

% Initialise upper bound based on x=zeros
x = zeros(2*n,1);
curOptimal = inf;

depth=1;

% Pop off the stack until it's empty
while (~isempty(STACK))
    
    % Pull off the first problem 
    curProblem = STACK(:,1);
    STACK(:,1) = [];
    
    % Calculate the optimisation to this problem
    [x_s,~,flag] = linprog(curProblem);
    
    % if it cannot be pruned by infeasibility or bound (i.e. is lower than
    % the current best legitimate candidate)
    if (~isempty(x_s) && c*x_s < curOptimal && flag == 1)
        
        % Need to create cell array of inventory vs injection for each
        % month
        d_s = x_s(1:end/2);
        e_s = x_s(end/2+1:end);
        
        % Calculate the inventory levels and the change in inventorylevels
        datapoints = [(V0+cumsum(e_s-d_s)*g) (e_s-d_s)*g];
        datapoints(:,1) = datapoints(:,1) - datapoints(:,2);
        datapoints = mat2cell(datapoints,ones(length(datapoints),1),2);
        
        relevantConstraints = {};
        
        for monthIndex=1:n
           relevantConstraints{monthIndex} = ...
               piecewiseConstraints{1 + (datapoints{monthIndex}(2) < 0)}(:,:,monthIndex);
           datapoints{monthIndex}(2) = abs(datapoints{monthIndex}(2)); 
        end 
        
        % Check against piecewise constraints
        [valid, invalidMonthIndex] = checkConstraints(datapoints, relevantConstraints);
       
        % If it satisfied the constraints (and is greater from before)
        if(valid)
            curOptimal = c*x_s;
            x = x_s;
            
        % It didn't satisfy constraints and is still greater, branch
        else
            
            % Subdivide the initial problem into two on either side of the point based
            % on the constraint that was violated (I or W, then which
            % segment)
            depth = depth+1;
            splitPoint = datapoints{invalidMonthIndex};
            
            subProblems = formSubproblems(  start, finish, curProblem, ...
                                            splitPoint, ...
                                            invalidMonthIndex);
            
            % breadth first
            % Reform the convex hull of the constraints for each subproblem
            % and add to the list
            % Need to use the splitpoint to re-relax
            [lowerConstraints, upperConstraints] = splitPiecewise(piecewiseConstraints, splitPoint);
            
            lowerProb = reformPiecewise(start, finish, cap, V0, subProblems{1}, lowerConstraints);
            upperProb = reformPiecewise(start, finish, cap, V0, subProblems{2}, upperConstraints);
            STACK = [STACK lowerProb upperProb];
        end
    end
end

fval = -curOptimal;

x(x<eps) = 0;
d = x(1:end/2);
e = x(end/2+1:end);
plotVariableConstraints(d,e,piecewiseConstraints,V0,Vn,L,U,cap)
return