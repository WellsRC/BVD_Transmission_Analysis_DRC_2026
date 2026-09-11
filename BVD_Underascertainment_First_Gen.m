clear;
clc;

N_Samp = 10^4;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);

Prob_Underascertainment_Samp = zeros(21,N_Samp);
avg_mu_Samp = zeros(21,N_Samp);

Prob_Underascertainment_F_Samp = zeros(21,N_Samp);
avg_mu_F_Samp = zeros(21,N_Samp);

Prob_Underascertainment_NF_Samp = zeros(21,N_Samp);
avg_mu_NF_Samp = zeros(21,N_Samp);

Prob_Underascertainment_temp = zeros(1,N_Samp);
avg_mu_temp = zeros(1,N_Samp);

Prob_Underascertainment_F_temp = zeros(1,N_Samp);
avg_mu_F_temp = zeros(1,N_Samp);

Prob_Underascertainment_NF_temp = zeros(1,N_Samp);
avg_mu_NF_temp = zeros(1,N_Samp);

for C1 = 0:20
    x1 = C1+[1:10^4];
    v_Underascertained = zeros(size(x1));
    v_avg = zeros(size(x1));
    parfor xx=1:length(x1)
        v_Underascertained(xx)=(1-integral(@(prob_asc)(betainc(prob_asc, C1 + 1, x1(xx) - C1)).*betapdf(prob_asc,max(C1,10^(-6)),x1(xx)-max(C1,10^(-6))),0,1));
        v_avg(xx) = 1-(1./v_Underascertained(xx)).*integral(@(prob_asc)(prob_asc-betainc(prob_asc, C1 + 1, x1(xx) - C1).*prob_asc).*betapdf(prob_asc,max(C1,10^(-6)),x1(xx)-max(C1,10^(-6))),0,1);
    end
    parfor ss=1:N_Samp
        % Baseline
        L_Underascertained=sum(nbinpdf(x1,k(ss),p(ss)).*v_Underascertained);
        avg_mu_temp(ss)=sum(nbinpdf(x1,k(ss),p(ss)).*v_avg)./sum(nbinpdf(x1,k(ss),p(ss)));
        L_Observed = nbincdf(C1,k(ss),p(ss));

        Prob_Underascertainment_temp(ss)=L_Underascertained./(L_Underascertained+L_Observed);

        % Fatal
        L_Underascertained=sum(nbinpdf(x1,k(ss),p_fatal(ss)).*v_Underascertained);
        avg_mu_F_temp(ss)=sum(nbinpdf(x1,k(ss),p_fatal(ss)).*v_avg)./sum(nbinpdf(x1,k(ss),p_fatal(ss)));
        L_Observed = nbincdf(C1,k(ss),p_fatal(ss));

        Prob_Underascertainment_F_temp(ss)=L_Underascertained./(L_Underascertained+L_Observed);


        % Non-fatal
        L_Underascertained=sum(nbinpdf(x1,k(ss),p_non_fatal(ss)).*v_Underascertained);
        avg_mu_NF_temp(ss)=sum(nbinpdf(x1,k(ss),p_non_fatal(ss)).*v_avg)./sum(nbinpdf(x1,k(ss),p_non_fatal(ss)));
        L_Observed = nbincdf(C1,k(ss),p_non_fatal(ss));

        Prob_Underascertainment_NF_temp(ss)=L_Underascertained./(L_Underascertained+L_Observed);
    end
    avg_mu_Samp (C1+1,:)=avg_mu_temp;
    avg_mu_F_Samp (C1+1,:)=avg_mu_F_temp;
    avg_mu_NF_Samp (C1+1,:)=avg_mu_NF_temp;

    Prob_Underascertainment_Samp(C1+1,:)=Prob_Underascertainment_temp;
    Prob_Underascertainment_F_Samp(C1+1,:)=Prob_Underascertainment_F_temp;
    Prob_Underascertainment_NF_Samp(C1+1,:)=Prob_Underascertainment_NF_temp;
end

Prob_Underascertainment= zeros(21,3);
avg_mu= ones(21,3);
avg_mu_Overall= ones(21,3);

Prob_Underascertainment_F= zeros(21,3);
avg_mu_F= ones(21,3);
avg_mu_F_Overall= ones(21,3);

Prob_Underascertainment_NF= zeros(21,3);
avg_mu_NF= ones(21,3);
avg_mu_NF_Overall= ones(21,3);

for gg=0:20
    beta_par=betafit(Prob_Underascertainment_Samp(gg+1,:));
    Prob_Underascertainment(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    Prob_Underascertainment(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
    
    beta_par=betafit(Prob_Underascertainment_F_Samp(gg+1,:));
    Prob_Underascertainment_F(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    Prob_Underascertainment_F(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));

    beta_par=betafit(Prob_Underascertainment_NF_Samp(gg+1,:));
    Prob_Underascertainment_NF(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    Prob_Underascertainment_NF(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
    
    if(gg>0)
        beta_par=betafit(avg_mu_Samp(gg+1,:));
        avg_mu(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
        avg_mu(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
        
        beta_par=betafit(avg_mu_F_Samp(gg+1,:));
        avg_mu_F(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
        avg_mu_F(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
    
        beta_par=betafit(avg_mu_NF_Samp(gg+1,:));
        avg_mu_NF(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
        avg_mu_NF(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
    end

    beta_par=betafit(avg_mu_Samp(gg+1,:).*Prob_Underascertainment_Samp(gg+1,:));
    avg_mu_Overall(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    avg_mu_Overall(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
    
    beta_par=betafit(avg_mu_F_Samp(gg+1,:).*Prob_Underascertainment_F_Samp(gg+1,:));
    avg_mu_F_Overall(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    avg_mu_F_Overall(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));

    beta_par=betafit(avg_mu_NF_Samp(gg+1,:).*Prob_Underascertainment_NF_Samp(gg+1,:));
    avg_mu_NF_Overall(gg+1,1)=(beta_par(1)-1)./(sum(beta_par)-2);
    avg_mu_NF_Overall(gg+1,2:3) = betainv([0.025 0.975],beta_par(1),beta_par(2));
end

save('Output_UA_1stG.mat','Prob_Underascertainment','avg_mu','avg_mu_Overall','Prob_Underascertainment_F','avg_mu_F','avg_mu_F_Overall','Prob_Underascertainment_NF','avg_mu_NF','avg_mu_NF_Overall')