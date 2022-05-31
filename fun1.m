function df=fun1(t,f,aplha,beta,e,C1,P_D,n,K)

PLR_min=C1(:,1); 
PLR_max=C1(:,2);
nP=n; ny=nP; 
nLam=1; nMu=ny+1; 
nInput=nP+ny+nLam+nMu; 
%beta=beta*0.1; 

% List the variables including introduced variables
p=f(1:nP); 
y=f(nP+1:nP+ny); 

%List the lagrangian multipliers for equality & inequality constraints
lambda=f(nP+ny+1:nP+ny+nLam); % Lagrangian inequality 
mu=f(nP+ny+nLam+1:nInput);  % Lagrangian equality 

% Calculate the diagonal matrixes
gamma_lambda=diag(lambda);gamma_mu=diag(mu.^2);

%dfp, dfy
dfp=3*C1(:,3).*(p.^2)+2*C1(:,4).*p+C1(:,5);
dfy=zeros(ny,1);


% consts 
g1=sum(y)-K;
g=[g1];
dg=[zeros(nP,1);
    ones(ny,1)];


h=[y.*(1-y);
    sum(p.*C1(:,7))-P_D]; 

dh=[zeros(nP,ny),C1(:,7);
    diag(1-2*y),zeros(ny,1)]; 


  
pyz=[p;y]-([dfp;dfy]+dg*lambda+dh*mu+aplha*dg*gamma_lambda*g+beta*dh*h);
py=min(1,max(0,pyz(nP+1:nP+ny))); 
pp=min(PLR_max.*py,max(PLR_min.*py,pyz(1:nP)));


df=[e*(-p+pp);
    e*(-y+py);
    e*(-lambda+max(0,lambda+g));
    e*(h)];
end