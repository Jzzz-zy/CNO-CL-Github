function [xn,vn]=pso(xp,v,N,pbest,gbest,nInput,PLR_max,n)

c1=0.5;%
c2=0.5;%
w=1; % Interia factor

p_max_check=[];
for j=1:N
    p_max_check=[p_max_check,PLR_max];
end

for j=1:N  
    r1=rand();
    r2=rand();
    vn(:,j)=w*v(:,j)+c1*r1*(pbest(1:nInput,j)-xp(:,j))+c2*r2*(gbest(1:nInput)-xp(:,j)); 
    xn(:,j)=xp(:,j)+vn(:,j); 
  
end
A=xn(1:n,:);
A(find(xn(1:n,:)>p_max_check))=p_max_check(find(xn(1:n,:)>p_max_check))*rand(1);
xn=[A;xn(n+1:nInput,:)];

end