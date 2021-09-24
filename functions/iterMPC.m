% Split future MPC
% ref: SOC target
% xk_1: present state
% mpcData: MPC struct

function [uk, mpcData] = iterMPC(ref,xk_1,mpcData)


    % Load MPC data
    Np = mpcData.Np;
    Nc = mpcData.Nc;
    uk_1 = mpcData.uk_1;
    Sigma = mpcData.Sigma;
    Ref = ref*ones(Np,1);
    Ru = mpcData.Ru;
    
%   Linearized State-space matrices
    Am = mpcData.matrices.A;
    Bm = mpcData.matrices.B;
    Csoc = mpcData.matrices.Cz;
    Dsoc = mpcData.matrices.Dz;
    

    % Augmented state vector
    dx = [xk_1;  uk_1];
    
   % Y = G*x + Phi*U
%  Compute SOC prediction matrices
   [augMatrices_soc]  = predMat(Am,Bm,Csoc,Dsoc,Nc,Np);
   [E,F,DU] = make_E_F_DU(augMatrices_soc,Np,Nc,Ru,Ref,dx);
   
    
   [M,gamma] = constraints(dx,mpcData);
    
   % Check if constraints are violated

    if sum(M*DU - gamma > 0) > 0
%         [DU] = QPhild(E,F,M,gamma);
        [DU] = hildreth(E,F,M,gamma,[],50);
    end
    
    du = DU(1);
    uk = du + uk_1;
    
    mpcData.uk_1 = uk;
    mpcData.SOCk_1 = xk_1(1);
    mpcData.DUk_1 = DU;

end

