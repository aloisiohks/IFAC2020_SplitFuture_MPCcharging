function [v,x,OCV] = iterModel(x, ik, model,deltaT)
    % Load constants
    
temp = 25 ; % Core temperature   Tc[k-1]    

R = getParamESC('RParam',temp,model);
RC = getParamESC('RCParam',temp,model);    
eta = getParamESC('etaParam',temp,model);  
Q = getParamESC('QParam',temp,model); 
R0 = getParamESC('R0Param',temp,model);

zk = x(1);
ir = x(2);


OCV = OCVfromSOCtemp(zk,25,model);
v = OCV  - R*ir - R0*ik;

A_RC = exp(-deltaT/RC);

A = diag([1 A_RC]);

B = [-deltaT*eta/(3600*Q)  ; ...
    (1-A_RC) ];
    
x = A*x + B*ik;



end