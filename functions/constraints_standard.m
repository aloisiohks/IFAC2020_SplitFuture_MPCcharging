function [ M, gamma ] = constraints_standard( dx, mpcData)

    Np = mpcData.Np;
    Nc = mpcData.Nc;
    Sigma = mpcData.Sigma;
    model = mpcData.model;  
    
    Am = mpcData.matrices.A;
    Bm = mpcData.matrices.B;
 
    
    C_z = mpcData.matrices.Cz;
    D_z = mpcData.matrices.Dz;
    C_v = mpcData.matrices.Cv;
    D_v = mpcData.matrices.Dv;

    
    [G_v,Phi_v] = predMat_standard(Am,Bm,C_v,D_v,Nc,Np);
    [G_z,Phi_z] = predMat_standard(Am,Bm,C_z,D_z,Nc,Np);
    
    if strcmp(mpcData.flag,'discharge')
    
    M = [ 
           Sigma;         % i  < i_max
          -Phi_z          % z > z_min]
          -Phi_v          % v  > v_min

        ] ;
    
     OCV =  OCVfromSOCtemp(dx(1),25,model);
     
    gamma = [ 
              mpcData.const.u_max*ones(Nc,1) - mpcData.uk_1*ones(Nc,1);         % i  < i_max
             -mpcData.const.z_min*ones(Np ,1) + G_z*dx ;                                     % z > z_min  
              (-mpcData.const.v_min*ones(Np ,1) + G_v*dx  + OCV*ones(Np ,1)  ) ;   % v  > v_min

            ];
    end
    
    
    if strcmp(mpcData.flag,'charge')
    
    M = [ 
           Phi_v;         % v  < v_max
          -Sigma;         % i  > i_min
           Phi_z;          % SOC<SOC_max
        ] ;
    
     OCV =  OCVfromSOCtemp(dx(1),25,model);
    
    gamma = [ 
              (mpcData.const.v_max*ones(Np ,1) - G_v*dx - OCV*ones(Np ,1)   );  % v  < v_max
             -mpcData.const.u_min*ones(Nc,1) + mpcData.uk_1*ones(Nc,1);         % i  > i_min
%               mpcData.const.z_max*ones(Np ,1) - G_z*dx - spkfData.zkbnd*ones(Np ,1);  % SOC<SOC_max
                         mpcData.const.z_max*ones(Np ,1) - G_z*dx ;  % SOC<SOC_max

              ];
    end
end
