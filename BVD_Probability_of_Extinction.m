clear;
clc;
N_Samp=10^6;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extinction at Generations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = [0:10];

Phi_V = zeros(length(G),N_Samp);

for ii=1:length(G)
    if(ii==1)
        Phi_V(ii,:)=p.^k;
    else
        Phi_V(ii,:)=p.^k./(1-Phi_V(ii-1,:)'.*(1-p)).^k;
    end
end

beta_par = betafit(Phi_V(G==0,:));

Phi_0 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_Phi_0=betainv(0.025,beta_par(1),beta_par(2));
ub_Phi_0=betainv(0.975,beta_par(1),beta_par(2));

clc;
fprintf(['Probabibility of extinction at the index case: ' num2str(Phi_0,'%4.3f') ' (95%% UI: ' num2str(lb_Phi_0,'%4.3f') char(8211) num2str(ub_Phi_0,'%4.3f') ') \n' ])

beta_par = betafit(Phi_V(G==5,:));
Phi_5 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_Phi_5=betainv(0.025,beta_par(1),beta_par(2));
ub_Phi_5=betainv(0.975,beta_par(1),beta_par(2));


fprintf(['Probabibility of extinction within 5 generations: ' num2str(Phi_5,'%4.3f') ' (95%% UI: ' num2str(lb_Phi_5,'%4.3f') char(8211) num2str(ub_Phi_5,'%4.3f') ') \n' ])

beta_par = betafit(Phi_V(G==10,:));
Phi_10 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_Phi_10=betainv(0.025,beta_par(1),beta_par(2));
ub_Phi_10=betainv(0.975,beta_par(1),beta_par(2));


fprintf(['Probabibility of extinction within 10 generations: ' num2str(Phi_10,'%4.3f') ' (95%% UI: ' num2str(lb_Phi_10,'%4.3f') char(8211) num2str(ub_Phi_10,'%4.3f') ') \n' ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fatal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_par = betafit(p_fatal.^k);

Phi_0 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_Phi_0=betainv(0.025,beta_par(1),beta_par(2));
ub_Phi_0=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Probabibility of extinction at the index case (Fatal): ' num2str(Phi_0,'%4.3f') ' (95%% UI: ' num2str(lb_Phi_0,'%4.3f') char(8211) num2str(ub_Phi_0,'%4.3f') ') \n' ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-fatal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_par = betafit(p_non_fatal.^k);

Phi_0 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_Phi_0=betainv(0.025,beta_par(1),beta_par(2));
ub_Phi_0=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Probabibility of extinction at the index case (Non-fatal): ' num2str(Phi_0,'%4.3f') ' (95%% UI: ' num2str(lb_Phi_0,'%4.3f') char(8211) num2str(ub_Phi_0,'%4.3f') ') \n' ])
