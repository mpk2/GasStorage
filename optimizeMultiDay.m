% This function will determine the optimal number of forward contracts to 
% buy / sell over a given period to optimise profits given constraints on 
% gas storage
%
% Inputs:
%   
%   n:  the number of months over which to optimise the contracts
%
%   f:  an rxn matrix where each row i represents the forward curve quoted
%       for a specific day d_i 1<=i<=r and each column specifies the quote 
%       for the month m_j 1<=j<=n such that f_i,j = forward contract price 
%       for month m_j in the forward curve for day d_i
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
%
%   fval:
%       the optimised value  
function [d,e,fval] = optimizeMultiDay(n, F, I, W, q, p, c, V0, Vn, L)

format short

s =[];
fval = 0;

for i=1:size(F,1)
    f = F(i,1:n);
    [d,e,fval] = optimizeContracts(n, f, I, W, q, p, c, V0, Vn, L);
    
    % Calculate c for this
    P = p(f) + c;
    Q = q(f) + c;
    c = [f-g*Q -f-g*P];
    
    % Add this amount to profit so far
    s(i,:) = [d e] - sum(s(1:i-1,:));
    fval = fval + c*s(i,:)';
end

d = s(:,1:(end/2));
e = s(:,(end/2):end);

end