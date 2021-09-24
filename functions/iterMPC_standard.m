% Standard MPC 
% ref: SOC target
% xk_1: present state
% mpcData: MPC struct

function [uk, mpcData] = iterMPC_standard(ref,xk_1,mpcData)


    % Load MPC data
    Np = mpcData.Np;
    Nc = mpcData.Nc;
    uk_1 = mpcData.uk_1;
    Sigma = mpcData.Sigma;
    Ref = ref*ones(Np,1);
    Ru = mpcData.Ru;
    
%   State Space matrices
    Am = mpcData.matrices.A;
    Bm = mpcData.matrices.B;
    Csoc = mpcData.matrices.Cz;
    Dsoc = mpcData.matrices.Dz;
    

    % Augmented state vector
    dx = [xk_1;  uk_1];
    
   % Y = G*x + Phi*U
%  Compute SOC prediction matrices
   [G_soc,Phi_soc]  = predMat_standard(Am,Bm,Csoc,Dsoc,Nc,Np);
   
   E = 2*(Phi_soc'*Phi_soc + Ru);
   F = -2*Phi_soc'*(Ref - G_soc*dx );
   DU = -E\F;
   
   [M,gamma] = constraints_standard(dx,mpcData);
    
   % Check if constraints are violated

    if sum(M*DU - gamma > 0) > 0
        [DU] = hildreth(E,F,M,gamma,[],200);
    end
    
    du = DU(1);
    uk = du + uk_1;
    
    mpcData.uk_1 = uk;
    mpcData.SOCk_1 = xk_1(1);
    mpcData.DUk_1 = DU;

end

