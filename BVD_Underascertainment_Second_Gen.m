clear;
clc;

N_Samp = 10^4;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);

Prob_Underascertainment_2 = zeros(31,21);
Prob_Underascertainment_2=Prob_Underascertainment_2(:);
avg_mu_2 = zeros(31,21);
avg_mu_2=avg_mu_2(:);

C1v=[0:20];
C2v=[0:30];

[C1v,C2v]=meshgrid(C1v,C2v);

C1v=C1v(:);
C2v=C2v(:);


for ii = 1:length(C1v)

    C1=C1v(ii);
    C2 = C2v(ii);
    x1 = C1v(ii)+[1:10^4];
    x2=C2v(ii)+[1:10^5];
    
    v_Underascertained_1 = zeros(size(x1));
    v_avg_1 = zeros(size(x1));
    parfor xx=1:length(x1)
        v_Underascertained_1(xx)=(1-integral(@(prob_asc)(betainc(prob_asc, C1 + 1, x1(xx) - C1)).*betapdf(prob_asc,max(C1,10^(-6)),x1(xx)-max(C1,10^(-6))),0,1));
        v_avg_1(xx) = (1-integral(@(prob_asc)(prob_asc+betainc(prob_asc, C1 + 1, x1(xx) - C1)-betainc(prob_asc, C1 + 1, x1(xx) - C1).*prob_asc).*betapdf(prob_asc,max(C1,10^(-6)),x1(xx)-max(C1,10^(-6))),0,1));
    end

    v_Underascertained_2 = zeros(size(x2));
    v_avg_2 = zeros(size(x2));
    parfor xx=1:length(x2)
        v_Underascertained_2(xx)=(1-integral(@(prob_asc)(betainc(prob_asc, C2 + 1, x2(xx) - C2)).*betapdf(prob_asc,max(C2,10^(-6)),x2(xx)-max(C2,10^(-6))),0,1));
        v_avg_2(xx) = (1-integral(@(prob_asc)(prob_asc+betainc(prob_asc, C2 + 1, x2(xx) - C2)-betainc(prob_asc, C2 + 1, x2(xx) - C2).*prob_asc).*betapdf(prob_asc,max(C2,10^(-6)),x2(xx)-max(C2,10^(-6))),0,1));
    end

    K= 1- sum(nbinpdf(x1,k_mle,p_mle).*nbincdf(C2,x1.*k_mle,p_mle));
       
    temp_1=nbinpdf(x1(:),k_mle,p_mle).*v_Underascertained_1(:);
    
    x1_temp=[1:C1v(ii)];    
    temp_x2_1=nbinpdf(repmat(x1_temp(:),1,length(x2)),k_mle,p_mle).*nbinpdf(repmat(x2,length(x1_temp),1),k_mle.*repmat(x1_temp(:),1,length(x2)),p_mle).*repmat(v_Underascertained_2,length(x1_temp),1);
    avg_temp_1 = temp_x2_1.*(repmat(x2,length(x1_temp),1).*repmat(v_avg_2,length(x1_temp),1))./(C1+repmat(x2,length(x1_temp),1));

    temp_x2_1=sum(temp_x2_1(:));
    avg_temp_1=sum(avg_temp_1(:));

    temp_x2_2=nbinpdf(x2,k_mle,p_mle).*v_Underascertained_2;
    
    temp_x2_3=nbinpdf(x1(:),k_mle,p_mle).*v_Underascertained_1(:).*nbincdf(C2v(ii),k_mle.*x1(:),p_mle);
    avg_temp_3 = temp_x2_3.*(x1(:).*v_avg_1(:)./(x1(:)+C2));
    temp_x2_3=sum(temp_x2_3(:));
    avg_temp_3=sum(avg_temp_3(:));

    temp_x2_2=temp_1*temp_x2_2;   
    avg_temp_2 = temp_x2_2.*(repmat(x1(:).*v_avg_1(:),1,length(x2))+repmat(x2.*v_avg_2,length(x1),1))./(repmat(x1(:),1,length(x2))+repmat(x2,length(x1),1));

    temp_x2_2=sum(temp_x2_2(:));
    avg_temp_2=sum(avg_temp_2(:));

    L_Underascertained = temp_x2_1+temp_x2_2+temp_x2_3;
    if(C1v(ii)>0)
        L_Observed = nbinpdf(0,k_mle,p_mle)+sum(nbinpdf(x1_temp,k_mle,p_mle).*nbincdf(C2v(ii),x1_temp.*k_mle,p_mle));
    elseif(C1v(ii)==0)
        if(C2v(ii)==0)
            L_Observed = nbinpdf(0,k_mle,p_mle);
        else
            L_Observed=0;
        end
    end
    
    Prob_Underascertainment_2(ii)=L_Underascertained./(L_Underascertained+L_Observed);
    avg_mu_2(ii)=(avg_temp_1+avg_temp_2+avg_temp_3)./K;
end

Prob_Underascertainment_2=reshape(Prob_Underascertainment_2,31,21);
avg_mu_2=reshape(avg_mu_2,31,21);

save('Output_UA_2ndG.mat','Prob_Underascertainment_2','avg_mu_2');
