function [x,v,pbest]=pnn(iteration,x,v,pbest,gbest,C1,P_D,n,N,PLR_max,PLR_min,K)

e=1000; 
aplha=10; 
beta=10;  

tend=0.3;
step=0.001;
tspan=0:step:tend;

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

parfor j =1:N
    f0=x(:,j);
    [t,f]=ode15s(@(t,y) fun1(t,y,aplha,beta,e,C1,P_D,n,K),tspan,f0,opts);
    [NRow,NCol]=size(f);
    f=f(NRow,:);
    xp(:,j)=f';
end
[nInput,nInputY]=size(xp);

% for j=1:N
%     while abs(sum(xp(1:n,j).*C1(:,7))- P_D)>1e-5
%         f0=rand(nInput,1);
%         [t,f]=ode15s(@(t,y) fun1(t,y,aplha,beta,e,C1,P_D,n,K),tspan,f0,opts);
%         [NRow,NCol]=size(f);
%         f=f(NRow,:);
%         xp(:,j)=f';
%     end
% end


if iteration ==1
    for j=1:N
        while abs(sum(xp(1:n,j).*C1(:,7))- P_D)>1e-5
            f0=[rand(n,1).*PLR_max;rand(nInput-n,1)];
            [t,f]=ode15s(@(t,y) fun1(t,y,aplha,beta,e,C1,P_D,n,K),tspan,f0,opts);
            [NRow,NCol]=size(f);
            f=f(NRow,:);
            xp(:,j)=f';
        end
    end
else
    for j=1:N
        if abs(sum(xp(1:n,j).*C1(:,7))- P_D)>1e-5
            xp(:,j)=pbest(1:nInput,j,iteration);
        end
    end    
end


for j=1:N
    for i=1:n
          Pij(i,j)=C1(i,3)*xp(i,j).^3+C1(i,4)*xp(i,j).^2+C1(i,5)*xp(i,j)+round(xp(i+n,j),2)*C1(i,6);
%         Pij(i,j)=C1(i,3)*xp(i,j).^3+C1(i,4)*xp(i,j).^2+C1(i,5)*xp(i,j)+C1(i,6);
    end
    if j==N
        ob=sum(Pij);
    end
end

if iteration>1
    for j=1:N
        if ob(j)<pbest(nInput+1,j,iteration)
            pbest(:,j,iteration+1)=[xp(:,j);ob(j)];
        else
            pbest(:,j,iteration+1)=pbest(:,j,iteration);
        end
    end
else
    pbest(:,:,iteration+1)=[xp;ob];
end

[x,v]=pso(xp,v,N,pbest(:,:,iteration),gbest(:,iteration),nInput,PLR_max,n);

end