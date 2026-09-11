clear;
clc;
% if isempty(gcp('nocreate'))
%     parpool(6);
% end
load('BVD_Model_Fit.mat','Parameters');

% k=Parameters(:,1);
% p=Parameters(:,2);
% mu = Parameters(:,4);
g_k =gamfit(Parameters(:,1));
k=(g_k(1)-1).*g_k(2);
beta_p = betafit(Parameters(:,2));
p= (beta_p(1)-1)./(sum(beta_p)-2);


mu = Parameters(:,4);
par_dist = fitdist(mu,'beta');
mu = (par_dist.a-1)./(par_dist.a+par_dist.b-2);

Num_Gen=17;
% uv = [0.05 0.1:0.1:0.9 0.95];
theta_G=zeros(length(k),Num_Gen+1); 
sigma_G_1=zeros(length(k),Num_Gen+1); 
sigma_G_2=zeros(length(k),Num_Gen+1); 
sigma_G_3=zeros(length(k),Num_Gen+1); 
sigma_G_4=zeros(length(k),Num_Gen+1); 
sigma_G_5=zeros(length(k),Num_Gen+1); 


for nn=1:length(k)
    theta_G(nn,1)=nbinpdf(0,k(nn),p(nn)); % Probability extinct w/ index
    
    for gg=2:Num_Gen+1
        theta_G(nn,gg)=p(nn).^k(nn)./(1-sum(theta_G(nn,1:(gg-1))).*(1-p(nn))).^k(nn)-sum(theta_G(nn,1:(gg-1)));
    end
end
for nn=1:length(k)
    sigma_G_1(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,mu(nn),k(nn),p(nn));
    sigma_G_2(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,2,mu(nn),k(nn),p(nn));
    sigma_G_3(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,3,mu(nn),k(nn),p(nn));
    sigma_G_4(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,4,mu(nn),k(nn),p(nn));
    sigma_G_5(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,5,mu(nn),k(nn),p(nn));
end

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

% Four  cases
pdf_4 =[sigma_G_4(:,1) diff(sigma_G_4,[],2)];
pdf_4=pdf_4./repmat(sum(pdf_4,2),1,Num_Gen+1);

pdf_4_trim = pdf_4(:,2:end);

E_4 = sum(theta_G.*pdf_4,2);

temp_EBD_4 = cumsum(theta_G,2).*(1-sigma_G_4);
temp_EBD_4=temp_EBD_4(:,1:end-1);
EBD_4 = sum(temp_EBD_4.*pdf_4_trim,2);

norm_pdf=temp_EBD_4.*pdf_4_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_4 =sum(Cases.*norm_pdf,2);

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