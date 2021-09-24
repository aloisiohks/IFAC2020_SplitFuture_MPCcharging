function [E,F,DU] = make_E_F_DU(augMatrices,Np,Nc,Ru,Ref,dx)

G = augMatrices.G;
Phi = augMatrices.Phi;
Gnf = augMatrices.Gnf;
Gff = augMatrices.Gff;
rw = Ru;

BarRs = ones(Np,1);
v = ones(1,Nc);


vT_GffT_Gff_v = v'*( Gff'* Gff)*v;
vT_GffT_Gnf = v'* Gff'* Gnf;
vT_GffT_R = v'* Gff'* BarRs ;
vT_GffT_Phi_At = v'* Gff'*G;
vT_GffT_Gff = v'* Gff'* Gff;
GnfT_Gnf = Gnf'* Gnf;
GnfT_Gff_v = Gnf'* Gff*v;
GnfT_Phi_At = Gnf'*G;
GnfT_R = Gnf'* BarRs ;
GnfT_Gff = Gnf'* Gff;

E = 2*( GnfT_Gnf + rw*eye(Nc ,Nc) - GnfT_Gff_v - vT_GffT_Gnf +...
vT_GffT_Gff_v );

F = -2*( GnfT_R *Ref(1) - vT_GffT_R *Ref(1) -...
GnfT_Phi_At*dx + vT_GffT_Phi_At * dx +...
GnfT_Gff*dx(end) - vT_GffT_Gff * dx(end) );

DU = inv( vT_GffT_Gff_v - vT_GffT_Gnf + GnfT_Gnf -...
GnfT_Gff_v + rw*eye(Nc ,Nc)) * ( GnfT_R *Ref(1) -...
vT_GffT_R *Ref(1) - GnfT_Phi_At * dx +...
vT_GffT_Phi_At * dx - vT_GffT_Gff * dx(end) +...
GnfT_Gff * dx(end) );

end