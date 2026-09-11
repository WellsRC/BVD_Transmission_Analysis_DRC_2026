clear;
clc;
% if isempty(gcp('nocreate'))
%     parpool(6);
% end
% 
% load("BVD_Model_Fit.mat",'Parameters')
% NS=size(Parameters,1);
% 
% R0=Parameters(:,1).*(1-Parameters(:,2))./Parameters(:,2);
% 
% k=Parameters(:,1);
% p=Parameters(:,2);
% u=Parameters(:,4);
% Num_Gen=15;
% 
% 
% theta_G=zeros(Num_Gen+1,length(k)); 
% sigma_G_1=zeros(Num_Gen+1,length(k)); 
% 
% parfor kk=1:length(k)
%     theta_G_temp=zeros(Num_Gen+1,1); 
%     theta_G_temp(1)=nbinpdf(0,k(kk),p(kk)); % Probability extinct w/ index
% 
%     for gg=2:Num_Gen+1
%         theta_G_temp(gg)=p(kk).^k(kk)./(1-sum(theta_G_temp(1:(gg-1))).*(1-p(kk))).^k(kk)-sum(theta_G_temp(1:(gg-1)));
%     end
%     theta_G(:,kk)=theta_G_temp;
%     sigma_G_1(:,kk)=Sample_Approximation_Detection_Chain(Num_Gen,1,u(kk),k(kk),p(kk));
% end
% 
% save('Output_Fig_S5.mat')

load('Output_Fig_S5.mat');
Gen_Det=zeros(1,length(k));
E_Det = zeros(1,length(k));
for kk=1:length(k)
    pdf_1 =[sigma_G_1(1,kk); diff(sigma_G_1(:,kk),[],1)];
    pdf_1 = pdf_1./sum(pdf_1);
    theta_G_1=theta_G(:,kk);
    Gen_Det(kk) = [0:Num_Gen]*pdf_1;
    E_Det(kk) = theta_G_1'*pdf_1;
end

OR = cumsum(theta_G,1).*(1-sigma_G_1)./(sigma_G_1.*(1-cumsum(theta_G,1)));

PRCC = zeros(3,size(OR,1));
for gen=0:Num_Gen
    X=[u OR(gen+1,:)'];
    Z=[k R0];
    [temp_r]=partialcorr(X,Z,'Type','Spearman');
    
    PRCC(1,gen+1)=temp_r(1,2);
    
    X=[R0 OR(gen+1,:)'];
    Z=[k u ];
    [temp_r]=partialcorr(X,Z,'Type','Spearman');
    
    PRCC(2,gen+1)=temp_r(1,2);
    
    
    X=[k OR(gen+1,:)'];
    Z=[u R0];
    [temp_r]=partialcorr(X,Z,'Type','Spearman');
    
    PRCC(3,gen+1)=temp_r(1,2);
end


close all;


figure('units','normalized','outerposition',[0.1 0.3 0.8 0.5])
subplot('Position',[0.075 0.15 0.9 0.8])
plot([0:Num_Gen],zeros(size([0:Num_Gen])),':','LineWidth',1.5,'Color',[0.5 0.5 0.5]); hold on;
p1=plot([0:Num_Gen],PRCC(1,:),'color',hex2rgb('#F5BE41'),'LineWidth',2,'Marker','o'); hold on;
p2=plot([0:Num_Gen],PRCC(2,:),'color',hex2rgb('#31A9B8'),'LineWidth',2,'Marker','square'); hold on;
p3=plot([0:Num_Gen],PRCC(3,:),'color',hex2rgb('#CF3721'),'LineWidth',2,'Marker','hexagram'); hold on;
box off;
set(gca,'tickdir','out','Fontsize',12,'LineWidth',2,'XTick',[0:Num_Gen],'YTick',[-1:0.2:1])
xlabel('Generation','FontSize',14)
ylabel('Partial Rank Correlation Coefficient','FontSize',14)
legend([p1 p2 p3],{'Underascertainment','Basic reproduction number','Dispersion coefficient'},'Position',[0.819907415427971,0.747004933051445,0.148123765635286,0.112050739957717],'FontSize',14)
theme('light');
print(gcf, 'Supplementary_Figure_5.png','-dpng','-r300');
