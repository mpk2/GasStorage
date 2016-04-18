function RealProblemSetup

    % Withdrawal rates (based on tight endpoint, e.g. 6400 mmbtu withdrawal
    % rate from 0 to 10% of inventory capacity)
    I = [   0       .35     .5      1; 
            9000    7500    6500    0];
        
    W = [   0       .1      .16     .35     1; 
            6400    8500    10300   12900   12900];    
    
    q = @(x) 0.015*x;
    p = @(x) 0;
    c = 0.01;
    
    % Calendar number for starting and finishing month (max 12 month)
    start = 4;
    finish = 3;
    
    % Doesn't need to be sent in
    N = 12;
    
    % Create random forward curve values
    F = 3+(25*[ceil(N/2):-1:1 1:floor(N/2)])/100;
    g=1e4;
    F=g*F';
   
    % Initial and final inventory levels
    V0 = 0;
    Vn = 0;
    
    % 1e6 mmbtu capacity facility
    cap=1000000;
    
    % Limit end of August to be 85% max (April-> August is 5th month)
    L = cap*[0 0 0 0 0 0 0 0 0 0 0 0];
    U = cap*[1 1 1 1 0.85 1 1 1 1 1 1 1];
    [d,e,fval]=optimizeContractsBB(start,finish,F,I,W,q,p,c,V0,Vn,L,U,cap)

end