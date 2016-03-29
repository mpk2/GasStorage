correct = 0;
f=[];
%for v=1:50
%for n=1:50

varlist = {'N','F','i','w','p','q','c','V0','Vn','l'};
clear(varlist{:});

N = 12;
F = 3+ceil(abs(5*randn(1,N)))/100;
F=F'*g;
i = abs(10000+500*randn(1,N));
w = abs(10000+500*randn(1,N));
p = @(x) 0.00;
q = @(x) 0.01;
c = 0;
V0 = 100000*randi(10)*rand(1);
Vn = ceil(2*rand*V0);
l = ceil(rand(1,N)*(V0+Vn)*0.2);

[d, e, fval] = optimizeContracts(1,12,F,i,w,q,p,c,V0,Vn,l);

%f(v,n) = fval;

if(~isempty([d' e' fval]))
    plotConstraints(d,e,i,w,V0,Vn,l);
end

%end
%end

% f = sum(f)/size(f,1);
% 
% figure;
% loglog(2.^((1:50)-1),f);
% title('Average Contribution of Injection and Withdrawal Constraints to Optimal Value');
% xlabel('Constraint Relaxation Factor')
% hold on;
% loglog(2.^((1:50)-1),(2.^((1:50)-1)).^(1))
% ylabel('Optimal Value of Objective')
% legend({'f', 'O(x)'},'Location','southeast')