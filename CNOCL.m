%% OCL
function [gbest,pbest,iteration,nInput,n]= CNOCL(P_D)
sysNo=4; 
if sysNo == 4
    load('FouChDa.mat');C1=FouChDa.C1;
    if P_D == 2610 || P_D == 2320 || P_D == 2030 || P_D == 1740
        K = 4;
        P=2;
        Max=10;
        cas=1;
    elseif P_D == 1450
        K=3;
        P=3;
        Max=10;
        cas=2;
    elseif P_D == 1160
        K=2;
        P=3;
        Max=10;
        cas=3;
    end
end


[n,nx]=size(C1); ep=0.1;

Div=0;
PLR_min=C1(:,1);
PLR_max=C1(:,2);

lower=PLR_min;
upper=PLR_max;

mut=0;mutup=0; mutlw=0;
values=FouChDa.Vs;
[nr,~]=size(values);
s=values(randi(nr),cas);rng(s);P=round(P*n/2+1);

times_end=1; iend=100;

opt_plot=[];nCount=1;

Po=P;
for ti=1:times_end
    iteration=1;
    gb_same=0;
    gb_thd=Max; 
    while iteration<=iend && gb_same<=gb_thd
        if iteration == 1 
            nP=n; ny=nP; 
            nLam=1; nMu=ny+1; 
            nInput=nP+ny+nLam+nMu;
            x=rand(nInput,Po);
            v=rand(nInput,Po);
            for j=1:Po
                for i=1:n
                    Pij(i,j)=C1(i,3)*x(i,j).^3+C1(i,4)*x(i,j).^2+C1(i,5)*x(i,j)+round(x(i+n,j),2)*C1(i,6); %x(i+n,j)*
                end
                if j==Po
                    Energy=sum(Pij);
                end
            end
            pbest=[x;Energy];
            [minObj,indx]=min(Energy);
            gbest(:,iteration)=pbest(:,indx);
        end

        [x,v,pbest]=pnn(iteration,x,v,pbest,gbest,C1,P_D,n,Po,PLR_max,PLR_min,K); % 75,70,65, k=5


        [mVal,ind_j]=min(pbest(nInput+1,:,iteration+1));
        if iteration > 1
            if mVal < gbest(nInput+1,iteration)
                gbest(:,iteration+1)=pbest(:,ind_j,iteration+1);
            else
                gbest(:,iteration+1)=gbest(:,iteration);
            end
        else
            gbest(:,iteration+1)=pbest(:,ind_j,iteration+1);
        end


        if abs(gbest(nInput+1,iteration+1)-gbest(nInput+1,iteration))<1e-4
            gb_same=gb_same+1;
        else
            gb_same=0;
        end

        Div=0;
        for j=1:Po
            Div_sin=norm(x(1:2*n,j)-gbest(1:2*n,iteration));
            Div=Div+Div_sin;
        end

        Div=Div/Po;

        if Div<ep
            a=exp(10*(iteration/iend));
            fai=-2.5*a+5*a*rand(1);
            eta=(1/sqrt(a))*exp(((-fai/a)^2)/2)*cos((5*(fai/a)));
            if eta>0
                for j=1:Po
                    x(1:n,j)=x(1:n,j)+eta*(upper-x(1:n,j));
                end
                iteration;
                mutup=mutup+1
            else
                for j=1:Po
                    x(1:n,j)=x(1:n,j)+eta*(x(1:n,j)-lower);
                end
                iteration;
                mutlw=mutlw+1
            end
        end
        iteration=iteration+1
    end

    mut=0;
    opt_plot=[opt_plot;Po,gbest(nInput+1,:)];
    bxplt=[opt_plot(:,1),opt_plot(:,iteration+1)];  % Save data for the boxplot

    iteration=iteration-1;

end

end



