% clear;
clc;

NS=5*10^3;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(NS);

p = p_fatal;
    uv=[0.01 0.1:0.1:0.9 0.99];

    Number_Unknown_Exposures=[1:10];
    PDF_test=zeros(length(Number_Unknown_Exposures),length(uv));
    lb_PDF_test=zeros(length(Number_Unknown_Exposures),length(uv));
    ub_PDF_test=zeros(length(Number_Unknown_Exposures),length(uv));
for uu=1:length(uv)
    u = uv(uu).*ones(NS,1);
    % Ratio of the the likelihood of human-to-human transmission relative to a
    % spillover event
    R=zeros(NS,1);
    
    current_G=1;
    x_past=[];
    for ss=1:NS
        zeta_G=zeros(10,1);
        for G=1:10
            if(G==1)
                zeta_G(G)=u(ss)-u(ss).*p(ss).^k(ss)./(1-u(ss).*(1-p(ss))).^k(ss);
            else
                temp_u = u(ss)-sum(zeta_G(1:G-1));
                zeta_G(G) = u(ss) -sum(zeta_G(1:G-1)) - u(ss).*p(ss).^k(ss)./(1-temp_u.*(1-p(ss))).^k(ss);
            end
        end
    
        L_HH=sum(zeta_G);
        L_S=1-u(ss);
    
        R(ss)=L_HH./(L_HH+L_S);
    end
    
    x_eval=linspace(0,1,101);
    for nn=1:length(Number_Unknown_Exposures)
        beta_HH_alt=betafit(1-(1-R).^nn);
        PDF_test(nn,uu)=(beta_HH_alt(1)-1)./(sum(beta_HH_alt)-2);
        lb_PDF_test(nn,uu)=betainv(0.025,beta_HH_alt(1),beta_HH_alt(2));
        ub_PDF_test(nn,uu)=betainv(0.975,beta_HH_alt(1),beta_HH_alt(2));
    end
end
