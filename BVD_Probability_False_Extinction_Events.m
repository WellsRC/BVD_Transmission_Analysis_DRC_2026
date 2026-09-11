clear;
clc;
close all;

N_Samp=10^4;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);
clc;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extinction vs under-ascertainment 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prob_Extinction_Index=p.^k;
Prob_Extinction_Index_F=p_fatal.^k;
Prob_Extinction_Index_NF=p_non_fatal.^k;

u =0.3;

Prob_Unascertained_First=p.^k.*(1-(1-u.*(1-p)).^k)./((1-u.*(1-p)).^k); 
beta_par = betafit(Prob_Unascertained_First);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

beta_par = betafit(Prob_Unascertained_First+Prob_Extinction_Index);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


P_samp = (1-u.*(1-p)).^k;
beta_par = betafit(P_samp);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

u =0.9;
Prob_Unascertained_First=p.^k.*(1-(1-u.*(1-p)).^k)./((1-u.*(1-p)).^k); 
beta_par = betafit(Prob_Unascertained_First);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

beta_par = betafit(Prob_Unascertained_First+Prob_Extinction_Index);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


P_samp = (1-u.*(1-p)).^k;

beta_par = betafit(P_samp);

P_90 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_90=betainv(0.025,beta_par(1),beta_par(2));
ub_P_90=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_90,'%4.3f') ' (95%% CrI: ' num2str(lb_P_90,'%4.3f') char(8211) num2str(ub_P_90,'%4.3f') ') \n \n' ])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extinction vs under-ascertainment: Fatal vs non-fatal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gp=gamfit(R0_Fatal);
R0_Fatal_mle = (gp(1)-1).*gp(2);
R0_Fatal_lb = gaminv(0.025,gp(1),gp(2));
R0_Fatal_ub = gaminv(0.975,gp(1),gp(2));
fprintf(['Fatal: Basic reporiduction number : ' num2str(R0_Fatal_mle,'%3.2f') ' (95%% CrI: ' num2str(R0_Fatal_lb,'%3.2f') char(8211) num2str(R0_Fatal_ub,'%3.2f') ') \n' ])

u =0.3;
Prob_Unascertained_First_F=p_fatal.^k.*(1-(1-u.*(1-p_fatal)).^k)./((1-u.*(1-p_fatal)).^k); 
beta_par = betafit(Prob_Unascertained_First_F);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Fatal: Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

beta_par = betafit(Prob_Unascertained_First_F+Prob_Extinction_Index_F);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Fatal: Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


P_samp = (1-u.*(1-p_fatal)).^k;

beta_par = betafit(P_samp);

P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));


fprintf(['Fatal: Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

u =0.9;
Prob_Unascertained_First_F=p_fatal.^k.*(1-(1-u.*(1-p_fatal)).^k)./((1-u.*(1-p_fatal)).^k); 
beta_par = betafit(Prob_Unascertained_First_F);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Fatal: Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


beta_par = betafit(Prob_Unascertained_First_F+Prob_Extinction_Index_F);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Fatal: Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


P_samp = (1-u.*(1-p_fatal)).^k;

beta_par = betafit(P_samp);

P_90 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_90=betainv(0.025,beta_par(1),beta_par(2));
ub_P_90=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Fatal: Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_90,'%4.3f') ' (95%% CrI: ' num2str(lb_P_90,'%4.3f') char(8211) num2str(ub_P_90,'%4.3f') ') \n \n' ])


gp=gamfit(R0_Non_fatal);
R0_Non_fatal_mle = (gp(1)-1).*gp(2);
R0_Non_fatal_lb = gaminv(0.025,gp(1),gp(2));
R0_Non_fatal_ub = gaminv(0.975,gp(1),gp(2));
fprintf(['Non-fatal: Basic reporiduction number : ' num2str(R0_Non_fatal_mle,'%3.2f') ' (95%% CrI: ' num2str(R0_Non_fatal_lb,'%3.2f') char(8211) num2str(R0_Non_fatal_ub,'%3.2f') ') \n' ])


u =0.3;
Prob_Unascertained_First_NF=p_non_fatal.^k.*(1-(1-u.*(1-p_non_fatal)).^k)./((1-u.*(1-p_non_fatal)).^k); 
beta_par = betafit(Prob_Unascertained_First_NF);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Non-fatal: Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


beta_par = betafit(Prob_Unascertained_First_NF+Prob_Extinction_Index_NF);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Non-fatal: Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])



P_samp = (1-u.*(1-p_non_fatal)).^k;

beta_par = betafit(P_samp);

P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));


fprintf(['Non-fatal: Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])

u =0.9;
Prob_Unascertained_First_NF=p_non_fatal.^k.*(1-(1-u.*(1-p_non_fatal)).^k)./((1-u.*(1-p_non_fatal)).^k); 
beta_par = betafit(Prob_Unascertained_First_NF);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Non-fatal: Probabibility of false extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


beta_par = betafit(Prob_Unascertained_First_NF+Prob_Extinction_Index_NF);
P_30 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_30=betainv(0.025,beta_par(1),beta_par(2));
ub_P_30=betainv(0.975,beta_par(1),beta_par(2));
fprintf(['Non-fatal: Probabibility of observing extinction at the index case (' num2str(u*100) '%%) : ' num2str(P_30,'%4.3f') ' (95%% CrI: ' num2str(lb_P_30,'%4.3f') char(8211) num2str(ub_P_30,'%4.3f') ') \n' ])


P_samp = (1-u.*(1-p_non_fatal)).^k;

beta_par = betafit(P_samp);

P_90 = (beta_par(1)-1)./(sum(beta_par)-2);

lb_P_90=betainv(0.025,beta_par(1),beta_par(2));
ub_P_90=betainv(0.975,beta_par(1),beta_par(2));

fprintf(['Non-fatal: Probabibility of extinction at the index case vs false extinction (' num2str(u*100) '%%) : ' num2str(P_90,'%4.3f') ' (95%% CrI: ' num2str(lb_P_90,'%4.3f') char(8211) num2str(ub_P_90,'%4.3f') ') \n' ])
