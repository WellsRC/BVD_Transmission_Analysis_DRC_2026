function [k,p,cfr,mu_A,I0,Generation,q_Parameter_MAT]=Sample_BVD_Parameters(k_OLD,p_OLD,cfr_OLD,mu_A_OLD,I0_OLD,Generation_OLD,weight)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Reproduction Number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gp_R0 = lsqnonlin(@(x) [1.46 1.68 1.36 1.83 1.24 2.15]-gaminv([0.35 0.65 0.2 0.8 0.05 0.95],x(1),x(2)),[2 2],[0 0],[100 100]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case fatality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://epiforecasts.io/BVDOutbreakSize/dev/analysis#Results (Accessed  July 23 2026)
bp_cfr = lsqnonlin(@(x) [0.44 0.49 0.41 0.53 0.36 0.59]-betainv([0.35 0.65 0.2 0.8 0.05 0.95],x(1),x(2)),[2 2],[0 0],[100 100]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dispersion coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_mle =0.27;
gp_k = lsqnonlin(@(x) [0.2 0.33] -gaminv([0.025 0.975],x(1),k_mle./(x(1)-1)),3,1,100);


if(sum(isnan(k_OLD))>0)
    q_Parameter_MAT = ones(6,1);
    
    R0=gamrnd(gp_R0(1),gp_R0(2),1,1);
   
    cfr = betarnd(bp_cfr(1),bp_cfr(2),1,1);    
    
    k = gamrnd(gp_k(1),k_mle./(gp_k(1)-1),1,1);
    
    p = k./(k+R0);

    mu_A = rand(1,1);
    I0 = randi(5,1,1);
    Generation = randi([1 17],1,1);
else
    % https://pmc.ncbi.nlm.nih.gov/articles/PMC2658655/
    dx=0.025;
    OLD_Parameter_MAT=[k_OLD,p_OLD,cfr_OLD,mu_A_OLD,I0_OLD,Generation_OLD];
    f_indx = find(rand(1)<=cumsum(weight),1,"first");
    OLD_Parameter_MAT = OLD_Parameter_MAT(f_indx,:);
    
    lb_Parameter_MAT = [OLD_Parameter_MAT(1:4)-dx OLD_Parameter_MAT(5:6)-1];   
    lb_Parameter_MAT([1:4])=max(lb_Parameter_MAT([1:4]),0);
    lb_Parameter_MAT([5:6])=max(lb_Parameter_MAT([5:6]),1);
    
    ub_Parameter_MAT = [OLD_Parameter_MAT(1:4)+dx OLD_Parameter_MAT(5:6)+1];   
    ub_Parameter_MAT([5])=min(ub_Parameter_MAT([5]),5);
    ub_Parameter_MAT([6])=min(ub_Parameter_MAT([6]),17);
    ub_Parameter_MAT([2:4])=min(ub_Parameter_MAT([2:4]),1);
    
    q_Parameter_MAT = 1./(ub_Parameter_MAT-lb_Parameter_MAT);
    
    
    Parameter_MAT = lb_Parameter_MAT+(ub_Parameter_MAT-lb_Parameter_MAT).*rand(size(ub_Parameter_MAT));
    k=Parameter_MAT(1);
    p=Parameter_MAT(2);
    cfr=Parameter_MAT(3);
    mu_A=Parameter_MAT(4);
    I0=round(Parameter_MAT(5));
    Generation=round(Parameter_MAT(6));
end
end