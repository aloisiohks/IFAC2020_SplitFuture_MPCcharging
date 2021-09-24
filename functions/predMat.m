function [augMatrices] = predMat(Am,Bm,Cm,D,Nc,Np)
% tic
% n = number of rows
% m = number of columms;
deltaN = 1;

[nA,mA] = size(Am);
[nB,mB] = size(Bm);
[nC,mC] = size(Cm);
[nD,mD] = size(D);

A = zeros(nA+mB,mA+mB);
A(1:nA,1:mA) = Am;
A(1:nA,nA+mB) = Bm;
A(nA+1:end,mA+1:end)=eye(mB);
B = zeros(nA+mB,1);
B(nA+mB,end) = 1;

C = zeros(1,mA+mB);
C(1,1:mC) = Cm;
C(1,mA+mB) = D;

G = zeros(Np,mA+mB);

for k=1:Np
    G(k,:) = C*A^(k);
    G2(k,:) = G(k,:)*A^(-1);
end


% for k=1:Np
%     G2(k,:) = C*A^(k-1);
% end

Phi = zeros(Np*nC,Nc*mB+deltaN);
for k=2:Nc+1+deltaN
    Phi(:,(k-2)+1:k-1)= [ zeros(k-2,1);G2(1:Np-(k-2),:)*B];
end

Gnf = Phi(:,1:Nc*mB);
Gff = Phi(:,Nc*mB+1:end);

augMatrices.A = A;
augMatrices.B = B;
augMatrices.C = C;
augMatrices.G = G;
augMatrices.Phi = Phi;
augMatrices.Gnf = Gnf;
augMatrices.Gff = Gff;

% toc
% disp('predMat')
end
