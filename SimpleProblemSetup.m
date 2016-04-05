function SimpleProblemSetup 

q = @(x) 0;
p = @(x) 0;
c = 0.01;
start = 3;
finish = 2;
N = 12;
F = 3+ceil(abs(5*randn(1,N)))/100;
g=1e4;
F=g*F';
I = [0 1; 9000 9000];
W = [0 1; 6400 6400];
V0 = 100;
Vn = 100;
cap=1000000;
L = cap*zeros(1,N);
[d,e,fval]=optimizeContractsBB(start,finish,F,I,W,q,p,c,V0,Vn,L,cap)
[d,e,fval]= optimizeContracts(start,finish, F, ...
    I(2,1)*ones(1,N), W(2,1)*ones(1,N), q, p, c, V0, Vn, L)

end