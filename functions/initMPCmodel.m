function [mpcData] = initMPCmodel(xhat,mpcData)

model = mpcData.model;
deltaT = mpcData.deltaT;

zkInd = 1; irkInd = 2;  

temp =25;
Q = getParamESC('QParam',temp,model); 
eta = getParamESC('etaParam',temp,model);


    
    R = getParamESC('RParam',temp,model);
    RC = getParamESC('RCParam',temp,model);    
    R0 = getParamESC('R0Param',temp,model); 
   
    Am = zeros(length(xhat),length(xhat));
    A_RC = exp(-deltaT/RC);
    Am(zkInd,zkInd) = 1;
    Am(irkInd,irkInd) = A_RC;


    B = zeros(length(xhat),1);
    B(zkInd) = -deltaT*eta/(3600*Q);
    B(irkInd) = 1 - A_RC;

    C_v = [0, -R];
    C_z = [1  0 ];

    Dv = -R0;
    Dz =0;


    
    mpcData.matrices.A = Am;
    mpcData.matrices.B = B;
    mpcData.matrices.Cv = C_v;
    mpcData.matrices.Cz = C_z;
    mpcData.matrices.Dv = Dv;
    mpcData.matrices.Dz = Dz;

end
