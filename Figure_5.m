clear;
clc;
% if isempty(gcp('nocreate'))
%     parpool(12);
% end
% load('BVD_Model_Fit.mat','Parameters');
% 
% k=Parameters(:,1);
% p=Parameters(:,2);
% mu = Parameters(:,4);
% cfr = Parameters(:,3);
% % g_k =gamfit(Parameters(:,1));
% % k_mle=(g_k(1)-1).*g_k(2);
% % beta_p = betafit(Parameters(:,2));
% % p_mle = (beta_p(1)-1)./(sum(beta_p)-2);
% % 
% % 
% % mu = Parameters(:,4);
% % par_dist = fitdist(mu,'beta');
% % mu_mle = (par_dist.a-1)./(par_dist.a+par_dist.b-2);
% 
% Num_Gen=17;
% % uv = [0.05 0.1:0.1:0.9 0.95];
% theta_G=zeros(length(k),Num_Gen+1); 
% sigma_G_1=zeros(length(k),Num_Gen+1); 
% Death_sigma_G_1=zeros(length(k),Num_Gen+1); 
% sigma_G_2=zeros(length(k),Num_Gen+1); 
% sigma_G_3=zeros(length(k),Num_Gen+1); 
% sigma_G_5=zeros(length(k),Num_Gen+1); 
% sigma_G_10=zeros(length(k),Num_Gen+1); 
% 
% for nn=1:length(k)
%     theta_G(nn,1)=nbinpdf(0,k(nn),p(nn)); % Probability extinct w/ index
% 
%     for gg=2:Num_Gen+1
%         theta_G(nn,gg)=p(nn).^k(nn)./(1-sum(theta_G(nn,1:(gg-1))).*(1-p(nn))).^k(nn)-sum(theta_G(nn,1:(gg-1)));
%     end
% end
% parfor nn=1:length(k)
%     sigma_G_1(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,mu(nn),k(nn),p(nn),1);
%     Death_sigma_G_1(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-cfr(nn),k(nn),p(nn),1);
%     sigma_G_2(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,2,mu(nn),k(nn),p(nn),1);
%     sigma_G_3(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,3,mu(nn),k(nn),p(nn),1);
%     sigma_G_5(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,5,mu(nn),k(nn),p(nn),1);
%     sigma_G_10(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,10,mu(nn),k(nn),p(nn),1);
% end
% 
% save('Extinction_Detection_Analysis.mat')

load('Extinction_Detection_Analysis.mat')

% Single death
pdf_d1 =[Death_sigma_G_1(:,1) diff(Death_sigma_G_1,[],2)];
pdf_d1=pdf_d1./repmat(sum(pdf_d1,2),1,Num_Gen+1);

pdf_d1_trim = pdf_d1(:,2:end);

E_d1 = sum(theta_G.*pdf_d1,2);

temp_EBD_1 = cumsum(theta_G,2).*(1-Death_sigma_G_1);
temp_EBD_1=temp_EBD_1(:,1:end-1);
EBD_d1 = sum(temp_EBD_1.*pdf_d1_trim,2);

% Single case
pdf_1 =[sigma_G_1(:,1) diff(sigma_G_1,[],2)];
pdf_1=pdf_1./repmat(sum(pdf_1,2),1,Num_Gen+1);

pdf_1_trim = pdf_1(:,2:end);

E_1 = sum(theta_G.*pdf_1,2);

temp_EBD_1 = cumsum(theta_G,2).*(1-sigma_G_1);
temp_EBD_1=temp_EBD_1(:,1:end-1);
EBD_1 = sum(temp_EBD_1.*pdf_1_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1.*pdf_1_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1 =sum(Cases.*norm_pdf,2);

% Two cases
pdf_2 =[sigma_G_2(:,1) diff(sigma_G_2,[],2)];
pdf_2=pdf_2./repmat(sum(pdf_2,2),1,Num_Gen+1);

pdf_2_trim = pdf_2(:,2:end);

E_2 = sum(theta_G.*pdf_2,2);

temp_EBD_2 = cumsum(theta_G,2).*(1-sigma_G_2);
temp_EBD_2=temp_EBD_2(:,1:end-1);
EBD_2 = sum(temp_EBD_2.*pdf_2_trim,2);

norm_pdf=temp_EBD_2.*pdf_2_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_2 =sum(Cases.*norm_pdf,2);

% Three cases
pdf_3 =[sigma_G_3(:,1) diff(sigma_G_3,[],2)];
pdf_3=pdf_3./repmat(sum(pdf_3,2),1,Num_Gen+1);

pdf_3_trim = pdf_3(:,2:end);

E_3 = sum(theta_G.*pdf_3,2);

temp_EBD_3 = cumsum(theta_G,2).*(1-sigma_G_3);
temp_EBD_3=temp_EBD_3(:,1:end-1);
EBD_3 = sum(temp_EBD_3.*pdf_3_trim,2);

norm_pdf=temp_EBD_3.*pdf_3_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_3 =sum(Cases.*norm_pdf,2);

% Five cases
pdf_5 =[sigma_G_5(:,1) diff(sigma_G_5,[],2)];
pdf_5=pdf_5./repmat(sum(pdf_5,2),1,Num_Gen+1);

pdf_5_trim = pdf_5(:,2:end);

E_5 = sum(theta_G.*pdf_5,2);

temp_EBD_5 = cumsum(theta_G,2).*(1-sigma_G_5);
temp_EBD_5=temp_EBD_5(:,1:end-1);
EBD_5 = sum(temp_EBD_5.*pdf_5_trim,2);

norm_pdf=temp_EBD_5.*pdf_5_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_5 =sum(Cases.*norm_pdf,2);

% Ten cases
pdf_10 =[sigma_G_10(:,1) diff(sigma_G_10,[],2)];
pdf_10=pdf_10./repmat(sum(pdf_10,2),1,Num_Gen+1);

pdf_10_trim = pdf_10(:,2:end);

E_10 = sum(theta_G.*pdf_10,2);

temp_EBD_10 = cumsum(theta_G,2).*(1-sigma_G_10);
temp_EBD_10=temp_EBD_10(:,1:end-1);
EBD_10 = sum(temp_EBD_10.*pdf_10_trim,2);

norm_pdf=temp_EBD_10.*pdf_10_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_10 =sum(Cases.*norm_pdf,2);


close all;
figure('units','normalized','outerposition',[0.2 0.2 0.4 0.7])


mle_pdf_1 =zeros(11,1);
mle_pdf_2 =zeros(11,1);
mle_pdf_3 =zeros(11,1);
mle_pdf_5 =zeros(11,1);
mle_pdf_10 =zeros(11,1);
mle_pdf_d1 =zeros(11,1);
mle_pdf_d1_bnds =zeros(11,2);
for gg = 0:10
    bp = betafit(pdf_1(:,gg+1));
    mle_pdf_1(gg+1)= (bp(1)-1)./(sum(bp)-2);
    if(gg>0)
        bp = betafit(pdf_2(:,gg+1));
        mle_pdf_2(gg+1)= (bp(1)-1)./(sum(bp)-2);
        
        bp = betafit(pdf_3(:,gg+1));
        mle_pdf_3(gg+1)= (bp(1)-1)./(sum(bp)-2);
    
        bp = betafit(pdf_5(:,gg+1));
        mle_pdf_5(gg+1)= (bp(1)-1)./(sum(bp)-2);
    
        bp = betafit(pdf_10(:,gg+1));
        mle_pdf_10(gg+1)= (bp(1)-1)./(sum(bp)-2);
    end

    bp = betafit(pdf_d1(:,gg+1));
    mle_pdf_d1(gg+1)= (bp(1)-1)./(sum(bp)-2);
    mle_pdf_d1_bnds(gg+1,:)=betainv([0.025 0.975],bp(1),bp(2));
end
Delay_Detection = zeros(size(pdf_1,1),5.*10^3);
for jj=1:size(pdf_1,1)

    cdf_1 = cumsum(pdf_1(jj,:))./sum(pdf_1(jj,:));
    cdf_d1 = cumsum(pdf_d1(jj,:))./sum(pdf_d1(jj,:));
    for ii=1:5.*10^3
        r=rand(1);
        Delay_Detection(jj,ii)= find(r<=cdf_d1,1,"first")-find(r<=cdf_1,1,"first");
    end
end
subplot('Position',[0.1 0.08 0.87 0.4])
plot([0:6],sqrt(mle_pdf_1(1:7)),'color',hex2rgb('#A3A599'),'LineWidth',2,'Marker','o','MarkerFaceColor',hex2rgb('#A3A599')); hold on
plot([0:6],sqrt(mle_pdf_d1(1:7)),'color',hex2rgb('#A10115'),'LineWidth',2,'Marker','s','MarkerFaceColor',hex2rgb('#A10115'))

set(gca,'LineWidth',2,'TickDir','out','XTick',[0:6],'FontSize',12,'Ytick',sqrt([0 0.01 0.05 0.1:.1:1]),'YTickLabel',{num2str([0 0.01 0.05 0.1:.1:1]')})
box off;
xlabel('Generation of detection','FontSize',14)
ylabel('Probability','FontSize',14)
legend({'One case','One death'},'FontSize',14)
legend boxoff;
xlim([0 6])
text(-0.1,1,'B','FontSize',28,'Units','normalized')
grid on;

subplot('Position',[0.1 0.56 0.87 0.4])
plot([0:10],sqrt(mle_pdf_1),'color',hex2rgb('#A3A599'),'LineWidth',2,'Marker','o','MarkerFaceColor',hex2rgb('#A3A599')); hold on
plot([0:10],sqrt(mle_pdf_2),'color',hex2rgb('#CE5A57'),'LineWidth',2,'Marker','d','MarkerFaceColor',hex2rgb('#CE5A57')); hold on
plot([0:10],sqrt(mle_pdf_3),'color',hex2rgb('#78A5A3'),'LineWidth',2,'Marker','pentagram','MarkerFaceColor',hex2rgb('#78A5A3')); hold on
plot([0:10],sqrt(mle_pdf_5),'color',hex2rgb('#E1B16A'),'LineWidth',2,'Marker','hexagram','MarkerFaceColor',hex2rgb('#E1B16A')); hold on
plot([0:10],sqrt(mle_pdf_10),'color',hex2rgb('#444C5C'),'LineWidth',2,'Marker','^','MarkerFaceColor',hex2rgb('#444C5C')); hold on

set(gca,'LineWidth',2,'TickDir','out','XTick',[0:10],'FontSize',12,'Ytick',sqrt([0 0.01 0.05 0.1:.1:1]),'YTickLabel',{num2str([0 0.01 0.05 0.1:.1:1]')})
box off;
legend({'One case','Two cases','Three cases','Five cases','Ten cases'},'FontSize',14)
legend boxoff;
xlabel('Generation of detection','FontSize',14)
ylabel('Probability','FontSize',14)
grid on;
xlim([0 10])
text(-0.1,1,'A','FontSize',28,'Units','normalized')
theme('light')

print(gcf,['Figure_3.png'],'-dpng','-r600');
% 
% clc;
% 
% bp_fit = betafit(E_1);
% mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
% bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% 
% display(['(Single Case) Probability of extinction with generation of detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% bp_fit = betafit(EBD_1);
% if(bp_fit(1)>1 && bp_fit(2)>1)
%     mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
%     bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% else
%     z = log(EBD_1./(1-EBD_1));
%     [m_temp,s_temp] = normfit(z);
%     mle = m_temp;
%     bnds =norminv([0.025 0.975],m_temp,s_temp);
%     mle = 1./(1+exp(-mle));
%     bnds = 1./(1+exp(-bnds));
% end
% 
% display(['(Single Case) Average probability of extinction  before detection detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% g_fit=gamfit(Num_Case_Missed_1-1);
% mle = 1+(g_fit(1)-1).*g_fit(2);
% bnds = gaminv([0.025 0.975],g_fit(1),g_fit(2))+1;
% display(['(Single Case) Number of cases in missed chains: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')'])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% bp_fit = betafit(E_2);
% mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
% bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% 
% display(['(Two Cases) Probability of extinction with generation of detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% bp_fit = betafit(EBD_2);
% if(bp_fit(1)>1 && bp_fit(2)>1)
%     mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
%     bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% else
%     z = log(EBD_2./(1-EBD_2));
%     [m_temp,s_temp] = normfit(z);
%     mle = m_temp;
%     bnds =norminv([0.025 0.975],m_temp,s_temp);
%     mle = 1./(1+exp(-mle));
%     bnds = 1./(1+exp(-bnds));
% end
% 
% display(['(Two Cases) Average probability of extinction  before detection detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% g_fit=gamfit(Num_Case_Missed_2-1);
% mle = 1+(g_fit(1)-1).*g_fit(2);
% bnds = gaminv([0.025 0.975],g_fit(1),g_fit(2))+1;
% display(['(Two Cases) Number of cases in missed chains: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')'])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% bp_fit = betafit(E_3);
% mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
% bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% 
% display(['(Three Cases) Probability of extinction with generation of detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% bp_fit = betafit(EBD_3);
% if(bp_fit(1)>1 && bp_fit(2)>1)
%     mle = (bp_fit(1)-1)./(sum(bp_fit)-2);
%     bnds = betainv([0.025 0.975],bp_fit(1),bp_fit(2));
% else
%     z = log(EBD_3./(1-EBD_3));
%     [m_temp,s_temp] = normfit(z);
%     mle = m_temp;
%     bnds =norminv([0.025 0.975],m_temp,s_temp);
%     mle = 1./(1+exp(-mle));
%     bnds = 1./(1+exp(-bnds));
% end
% 
% display(['(Three Cases) Average probability of extinction  before detection detection: ' num2str(mle,'%4.3f') ' (95% CrI: ' num2str(bnds(1),'%4.3f') char(8211) num2str(bnds(2),'%4.3f') ')'])
% 
% g_fit=gamfit(Num_Case_Missed_3-1);
% mle = 1+(g_fit(1)-1).*g_fit(2);
% bnds = gaminv([0.025 0.975],g_fit(1),g_fit(2))+1;
% display(['(Three Cases) Number of cases in missed chains: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')'])