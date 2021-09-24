function [ M, gamma ] = constraints( dx, mpcData)

    Np = mpcData.Np;
    Nc = mpcData.Nc;
    Sigma = mpcData.Sigma;
    model = mpcData.model;  
    
    Am = mpcData.matrices.A;
    Bm = mpcData.matrices.B;
 
    
    Csoc = mpcData.matrices.Cz;
    Dsoc = mpcData.matrices.Dz;
    Cv = mpcData.matrices.Cv;
    Dv = mpcData.matrices.Dv;

    
    [augMatrices_v]  = predMat(Am,Bm,Cv,Dv,Nc,Np);
    [augMatrices_soc]  = predMat(Am,Bm,Csoc,Dsoc,Nc,Np);
    G_v = augMatrices_v.G;
    G_z = augMatrices_soc.G;
    Gnf_v = augMatrices_v.Gnf;
    Gnf_z = augMatrices_soc.Gnf;
    Gff_v = augMatrices_v.Gff;
    Gff_z = augMatrices_soc.Gff;
    Phi_v = augMatrices_v.Phi;
    Phi_z = augMatrices_soc.Phi;
    
    v = ones(1,Nc);
    
    if strcmp(mpcData.flag,'discharge')
    
    M = [ 
           Sigma;         % i  < i_max
          -( Gnf_z - Gff_z *v);         % z > z_min]
          -( Gnf_v - Gff_v *v);          % v  > v_min

        ] ;
    
     OCV =  OCVfromSOCtemp(dx(1),25,model);
     
    gamma = [ 
              mpcData.const.u_max*ones(Nc,1) - mpcData.uk_1*ones(Nc,1);         % i  < i_max
             -mpcData.const.z_min*ones(Np ,1) + G_z*dx -Gff_z*dx(end) ;                                     % z > z_min  
              (-mpcData.const.v_min*ones(Np ,1) + G_v*dx  + OCV*ones(Np ,1) - Gff_v*dx(end) ) ;   % v  > v_min

            ];
    end
    
    
    if strcmp(mpcData.flag,'charge')
    
    M = [ 
           ( Gnf_v - Gff_v*v);         % v  < v_max
          -Sigma;         % i  > i_min
           ( Gnf_z - Gff_z *v);          % SOC<SOC_max
        ] ;
    
     OCV =  OCVfromSOCtemp(dx(1),25,model);
    
    gamma = [ 
              (mpcData.const.v_max*ones(Np ,1) - G_v*dx - OCV*ones(Np ,1) + Gff_v*dx(end)  );  % v  < v_max
             -mpcData.const.u_min*ones(Nc,1) + mpcData.uk_1*ones(Nc,1);         % i  > i_min
%               mpcData.const.z_max*ones(Np ,1) - G_z*dx - spkfData.zkbnd*ones(Np ,1);  % SOC<SOC_max
                         mpcData.const.z_max*ones(Np ,1) - G_z*dx + Gff_z*dx(end) ;  % SOC<SOC_max

              ];
    end
end
