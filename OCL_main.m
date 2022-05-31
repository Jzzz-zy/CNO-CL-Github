%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo code of the paper <Optimal Chiller Loading Based on Collaborative 
% Neurodynamic Optimization>, 
% accepted by IEEE Transactions on Industrial Informatics. 
% 
% [Example]: Four-chiller system 
% 
% Contact CHEN Zhongying (zy.chen@my.cityu.edu.hk) if you have
% any questions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc
rng(1)
%% load
P_D_Pool=[2610, 2320, 2030 ,1740, 1450, 1160];
P_D=1160;

%% compu
[gbest,pbest,iteration,nInput,n]= CNOCL(P_D);


%% Plot
figure()
tiledlayout(1,2)

nexttile
plot(0:iteration,gbest(nInput+1,1:iteration+1))
xlabel('Iteration')
ylabel('Power consumption (kW)')

nexttile
if n==iteration+1
    gbest_re=gbest(1:n,:)';
else
    gbest_re=gbest(1:n,:);
end
plot(0:iteration,gbest_re(1:n,:))
xlabel('Iteration')
ylabel('PLR')


Power=gbest(end,end);
fprintf(' Load: %.0f RT; \n Power consumption: %.2f kW \n',P_D,Power)
%PLR=round(gbest(1:n,end)',2);
%display(PLR)
