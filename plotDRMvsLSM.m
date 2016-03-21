varlist = {'N','F','i','w','p','q','c','V0','Vn','l'};
clear(varlist{:});

N = 12;
p = @(x) 0.00;
q = @(x) 0.00;
c = 0;
V0 = 100000*randi(10)*rand(1);
Vn = ceil(rand*V0);
L = ceil(randi(40,1,N)*(V0+Vn)/200);
fvalSingle = [];
fvalMulti = [];

F=[];
i = abs(500*rand(1,N));
w = abs(500*rand(1,N));


for month=1:12
    
    F(end+1,:) = (3+ceil(abs(5*randn(1,N))))*1e2
   
    [d, e, fvalSingle(end+1)] = optimizeContracts(N, F(end,:),i,w,q,p,c,V0,Vn,L);
   
end

[d, e, fvalMulti] = optimizeMultiDay(N,F,i,w,q,p,c,V0,Vn,L);

fvalSingle
fvalMulti
plot(1:12, fvalSingle);
hold on
plot(1:12, fvalMulti(2:end));

