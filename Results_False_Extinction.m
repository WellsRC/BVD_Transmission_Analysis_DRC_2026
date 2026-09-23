clear;
clc;
close all;

load("BVD_Model_Fit.mat",'Parameters')
k=Parameters(:,1);
p=Parameters(:,2);
u=Parameters(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probability of extinction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = optimoptions("fmincon","Display","none","FunctionTolerance",10^(-9));
Prob_True_Extinction=p.^k;
Prob_False_Extinction=p.^k.*(1-(1-u.*(1-p)).^k)./((1-u.*(1-p)).^k); 
Prob_Apparent_Extinction = Prob_True_Extinction + Prob_False_Extinction;
Conditional_Prob_Extinction = (1-u.*(1-p)).^k;
% False extinction
b = fitdist(Prob_False_Extinction,'Kernel','Support',[0 1]);
mle = fmincon(@(z)-log(pdf(b,z)),mean(Prob_False_Extinction),[],[],[],[],0,1,[],opts);
mle_false_ext = mle;
lb = icdf(b,0.025);
ub = icdf(b,0.975);

display(['The probability of false extinction: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')']);

b = fitdist(1./Prob_False_Extinction,'Kernel','Support','positive');
mle = fmincon(@(z)-log(pdf(b,z)),mean(1./Prob_False_Extinction),[],[],[],[],0,100,[],opts);
lb = icdf(b,0.025);
ub = icdf(b,0.975);

display(['The average number of ascertained cases for false-extinction event to occurr: ' num2str(round(mle),'%d') ' (95% CrI: ' num2str(round(lb),'%d') char(8211) num2str(round(ub),'%d') ')']);

% Apparent extinction
b = fitdist(Prob_Apparent_Extinction,'Kernel','Support',[0 1]);
mle = fmincon(@(z)-log(pdf(b,z)),mean(Prob_Apparent_Extinction),[],[],[],[],0,1,[],opts);
lb = icdf(b,0.025);
ub = icdf(b,0.975);

display(['The probability of apparent extinction: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')']);

% True extinction: Conditional on extinction
b = fitdist(Conditional_Prob_Extinction,'Kernel','Support',[0 1]);
mle = fmincon(@(z)-log(pdf(b,z)),mean(Conditional_Prob_Extinction),[],[],[],[],0,1,[],opts);
lb = icdf(b,0.025);
ub = icdf(b,0.975);

display(['Conditional on apparent extinctio, the probability of true extinction: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')']);

N_Confirmed = 5514;
x=[0:N_Confirmed];

False_Ext_Events = zeros(length(k),1);

L = binopdf(x,N_Confirmed,mle_false_ext);
False_Ext_Events_mle = x(L==max(L));

for jj=1:length(False_Ext_Events)
    L = binopdf(x,N_Confirmed,Prob_False_Extinction(jj));
    False_Ext_Events(jj)=x(L==max(L));
end

lb = prctile(False_Ext_Events,2.5);
ub = prctile(False_Ext_Events,97.5);

display(['From ' num2str(N_Confirmed) ' cases, the number of false extinction events: ' num2str(False_Ext_Events_mle,'%d') ' (95% CrI: ' num2str(lb,'%d') char(8211) num2str(ub,'%d') ')']);

