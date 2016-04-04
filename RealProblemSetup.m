function RealProblemSetup

W = [0 .1 .16 .35; 6400 8500 10300 12900];
q = @(x) 0.00;
p = @(x) 0;
c = 0.01;
start = 4;
finish = 3;
N = 12;
F = 3+ceil(50*randn(1,N))/100;
g=1e4;
F=g*F';
I = [0 .35 .5 1; 9000 7500 6500 0];
W = [0 .1 .16 .35 1; 6400 8500 10300 12900 12900];
V0 = 1000;
Vn = 1000;
cap=1000000;
L = cap*[0 0 0 0 1 0 0 0 0 0 0 0];
[d,e,fval]=optimizeContractsBB(start,finish,F,I,W,q,p,c,V0,Vn,L,cap)

end