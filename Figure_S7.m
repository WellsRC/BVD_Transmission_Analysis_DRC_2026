% clear;
% clc;
% if isempty(gcp('nocreate'))
%     parpool(16);
% end
% 
% rng(20260720)
% 
% [~,~,~,k,p,cfr,~,~] = BVD_Parameters(5000);
% 
% Num_Gen=17;
% theta_G=zeros(length(k),Num_Gen+1); 
% 
% sigma_G_1_A10=zeros(length(k),Num_Gen+1); 
% sigma_G_1_A30=zeros(length(k),Num_Gen+1); 
% sigma_G_1_A50=zeros(length(k),Num_Gen+1);
% sigma_G_1_A70=zeros(length(k),Num_Gen+1); 
% sigma_G_1_A90=zeros(length(k),Num_Gen+1); 
% 
% for nn=1:length(k)
%     theta_G(nn,1)=nbinpdf(0,k(nn),p(nn)); % Probability extinct w/ index
% 
%     for gg=2:Num_Gen+1
%         theta_G(nn,gg)=p(nn).^k(nn)./(1-sum(theta_G(nn,1:(gg-1))).*(1-p(nn))).^k(nn)-sum(theta_G(nn,1:(gg-1)));
%     end
% end
% parfor nn=1:length(k)
%     sigma_G_1_A10(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-0.1,k(nn),p(nn),1);
%     sigma_G_1_A30(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-0.3,k(nn),p(nn),1);
%     sigma_G_1_A50(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-0.5,k(nn),p(nn),1);
%     sigma_G_1_A70(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-0.7,k(nn),p(nn),1);
%     sigma_G_1_A90(nn,:)=Sample_Approximation_Detection_Chain(Num_Gen,1,1-0.9,k(nn),p(nn),1);
% end
% 
% save('Extinction_Detection_Analysis_Figure_S7.mat')

load('Extinction_Detection_Analysis_Figure_S7.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ascertainment rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10% Ascertainemnt 
pdf_1_A10 =[sigma_G_1_A10(:,1) diff(sigma_G_1_A10,[],2)];
pdf_1_A10=pdf_1_A10./repmat(sum(pdf_1_A10,2),1,Num_Gen+1);

pdf_1_A10_trim = pdf_1_A10(:,2:end);

E_1_A10 = sum(theta_G.*pdf_1_A10,2);

temp_EBD_1_A10 = cumsum(theta_G,2).*(1-sigma_G_1_A10);
temp_EBD_1_A10=temp_EBD_1_A10(:,1:end-1);
EBD_1_A10 = sum(temp_EBD_1_A10.*pdf_1_A10_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1_A10.*pdf_1_A10_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1_A10 =sum(Cases.*norm_pdf,2);

% 30% Ascertainemnt 
pdf_1_A30 =[sigma_G_1_A30(:,1) diff(sigma_G_1_A30,[],2)];
pdf_1_A30=pdf_1_A30./repmat(sum(pdf_1_A30,2),1,Num_Gen+1);

pdf_1_A30_trim = pdf_1_A30(:,2:end);

E_1_A30 = sum(theta_G.*pdf_1_A30,2);

temp_EBD_1_A30 = cumsum(theta_G,2).*(1-sigma_G_1_A30);
temp_EBD_1_A30=temp_EBD_1_A30(:,1:end-1);
EBD_1_A30 = sum(temp_EBD_1_A30.*pdf_1_A30_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1_A30.*pdf_1_A30_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1_A30 =sum(Cases.*norm_pdf,2);

% 50% Ascertainemnt 
pdf_1_A50 =[sigma_G_1_A50(:,1) diff(sigma_G_1_A50,[],2)];
pdf_1_A50=pdf_1_A50./repmat(sum(pdf_1_A50,2),1,Num_Gen+1);

pdf_1_A50_trim = pdf_1_A50(:,2:end);

E_1_A50 = sum(theta_G.*pdf_1_A50,2);

temp_EBD_1_A50 = cumsum(theta_G,2).*(1-sigma_G_1_A50);
temp_EBD_1_A50=temp_EBD_1_A50(:,1:end-1);
EBD_1_A50 = sum(temp_EBD_1_A50.*pdf_1_A50_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1_A50.*pdf_1_A50_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1_A50 =sum(Cases.*norm_pdf,2);

% 70% Ascertainemnt 
pdf_1_A70 =[sigma_G_1_A70(:,1) diff(sigma_G_1_A70,[],2)];
pdf_1_A70=pdf_1_A70./repmat(sum(pdf_1_A70,2),1,Num_Gen+1);

pdf_1_A70_trim = pdf_1_A70(:,2:end);

E_1_A70 = sum(theta_G.*pdf_1_A70,2);

temp_EBD_1_A70 = cumsum(theta_G,2).*(1-sigma_G_1_A70);
temp_EBD_1_A70=temp_EBD_1_A70(:,1:end-1);
EBD_1_A70 = sum(temp_EBD_1_A70.*pdf_1_A70_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1_A70.*pdf_1_A70_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1_A70 =sum(Cases.*norm_pdf,2);

% 90% Ascertainemnt 
pdf_1_A90 =[sigma_G_1_A90(:,1) diff(sigma_G_1_A90,[],2)];
pdf_1_A90=pdf_1_A90./repmat(sum(pdf_1_A90,2),1,Num_Gen+1);

pdf_1_A90_trim = pdf_1_A90(:,2:end);

E_1_A90 = sum(theta_G.*pdf_1_A90,2);

temp_EBD_1_A90 = cumsum(theta_G,2).*(1-sigma_G_1_A90);
temp_EBD_1_A90=temp_EBD_1_A90(:,1:end-1);
EBD_1_A90 = sum(temp_EBD_1_A90.*pdf_1_A90_trim,2);

R0 = k.*(1-p)./p;
Cases =cumsum(repmat(R0,1,Num_Gen).^repmat([0:Num_Gen-1],length(k),1),2);
norm_pdf=temp_EBD_1_A90.*pdf_1_A90_trim;
norm_pdf = norm_pdf./repmat(sum(norm_pdf,2),1,Num_Gen);
Num_Case_Missed_1_A90 =sum(Cases.*norm_pdf,2);

close all;

opts = optimoptions("fmincon","Display","none","FunctionTolerance",10^(-9));


mle_pdf_1_A10 =zeros(11,1);
mle_pdf_1_A30 =zeros(11,1);
mle_pdf_1_A50 =zeros(11,1);
mle_pdf_1_A70 =zeros(11,1);
mle_pdf_1_A90 =zeros(11,1);
xx_temp =linspace(0,1,1001);
for gg = 0:10
   
    if(gg>0)
        bp = fitdist(pdf_1_A10(:,gg+1),'Kernel','Support',[-10^(-16) 1]);
        y=pdf(bp,xx_temp);
        if(max(y)>0 && max(y)~=y(1))
            mle_pdf_1_A10(gg+1)= fmincon(@(z)-log(pdf(bp,z)),mean(pdf_1_A10(:,gg+1)),[],[],[],[],0,1,[],opts);
        end

        bp = fitdist(pdf_1_A30(:,gg+1),'Kernel','Support',[-10^(-16) 1]);
        y=pdf(bp,xx_temp);
        if(max(y)>0 && max(y)~=y(1))
            mle_pdf_1_A30(gg+1)= fmincon(@(z)-log(pdf(bp,z)),mean(pdf_1_A30(:,gg+1)),[],[],[],[],0,1,[],opts);
        end
    

        bp = fitdist(pdf_1_A50(:,gg+1),'Kernel','Support',[-10^(-16) 1]);
        y=pdf(bp,xx_temp);
        if(max(y)>0 && max(y)~=y(1))
            mle_pdf_1_A50(gg+1)= fmincon(@(z)-log(pdf(bp,z)),mean(pdf_1_A50(:,gg+1)),[],[],[],[],0,1,[],opts);
        end
    
        bp = fitdist(pdf_1_A70(:,gg+1),'Kernel','Support',[-10^(-16) 1]);
        y=pdf(bp,xx_temp);
        if(max(y)>0 && max(y)~=y(1))
            mle_pdf_1_A70(gg+1)= fmincon(@(z)-log(pdf(bp,z)),mean(pdf_1_A70(:,gg+1)),[],[],[],[],0,1,[],opts);
        end
    
    
        bp = fitdist(pdf_1_A90(:,gg+1),'Kernel','Support',[-10^(-16) 1]);
        y=pdf(bp,xx_temp);
        if(max(y)>0 && max(y)~=y(1))
            mle_pdf_1_A90(gg+1)= fmincon(@(z)-log(pdf(bp,z)),mean(pdf_1_A90(:,gg+1)),[],[],[],[],0,1,[],opts);
        end

    else
        mle_pdf_1_A10(gg+1)= mean(pdf_1_A10(:,gg+1));
        mle_pdf_1_A30(gg+1)= mean(pdf_1_A30(:,gg+1));
        mle_pdf_1_A50(gg+1)= mean(pdf_1_A50(:,gg+1));  
        mle_pdf_1_A70(gg+1)= mean(pdf_1_A70(:,gg+1));
        mle_pdf_1_A90(gg+1)= mean(pdf_1_A90(:,gg+1));

    end

end
figure('units','normalized','outerposition',[0.2 0.2 0.4 0.45])
subplot('Position',[0.1 0.14 0.87 0.8])
cond_mle_pdf_1_A10 = mle_pdf_1_A10(2:end)./sum(mle_pdf_1_A10(2:end));
cond_mle_pdf_1_A30 = mle_pdf_1_A30(2:end)./sum(mle_pdf_1_A30(2:end));
cond_mle_pdf_1_A50 = mle_pdf_1_A50(2:end)./sum(mle_pdf_1_A50(2:end));
cond_mle_pdf_1_A70 = mle_pdf_1_A70(2:end)./sum(mle_pdf_1_A70(2:end));
cond_mle_pdf_1_A90 = mle_pdf_1_A90(2:end)./sum(mle_pdf_1_A90(2:end));

plot([1:10],sqrt(cond_mle_pdf_1_A10),'color',hex2rgb('#CE5A57'),'LineWidth',2,'Marker','d','MarkerFaceColor',hex2rgb('#CE5A57')); hold on
plot([1:10],sqrt(cond_mle_pdf_1_A30),'color',hex2rgb('#78A5A3'),'LineWidth',2,'Marker','pentagram','MarkerFaceColor',hex2rgb('#78A5A3')); hold on
plot([1:10],sqrt(cond_mle_pdf_1_A50),'color',hex2rgb('#E1B16A'),'LineWidth',2,'Marker','hexagram','MarkerFaceColor',hex2rgb('#E1B16A')); hold on
plot([1:10],sqrt(cond_mle_pdf_1_A70),'color',hex2rgb('#A3A599'),'LineWidth',2,'Marker','o','MarkerFaceColor',hex2rgb('#A3A599')); hold on
plot([1:10],sqrt(cond_mle_pdf_1_A90),'color',hex2rgb('#444C5C'),'LineWidth',2,'Marker','^','MarkerFaceColor',hex2rgb('#444C5C')); hold on

set(gca,'LineWidth',2,'TickDir','out','XTick',[0:10],'FontSize',12,'Ytick',sqrt([0 0.01 0.05 0.1:.1:1]),'YTickLabel',{num2str([0 0.01 0.05 0.1:.1:1]')})
box off;
legend({'10% ascertainment','30% ascertainment','50% ascertainment','70% ascertainment','90% ascertainment'},'FontSize',14)
legend boxoff;
xlabel('Generation of detection','FontSize',14)
ylabel('Probability','FontSize',14)
grid on;
xlim([1 10])
ylim([0 1])
theme('light')
print(gcf,['Figure_S7.png'],'-dpng','-r600');
clc;

% Generation: 10% Ascertainment 
Gv=[1:10];
G_mle = Gv(cond_mle_pdf_1_A10==max(cond_mle_pdf_1_A10));
G_95 = Gv(find(cumsum(cond_mle_pdf_1_A10)>0.95,1,"first"));

display(['With 10% ascertainment and conditional index not ascertained, the Generation of detection: ' num2str(G_mle,'%d') ' (95 percentile = ' num2str(G_95,'%d') ')']);

% Generation: 30% Ascertainment 
Gv=[1:10];
G_mle = Gv(cond_mle_pdf_1_A30==max(cond_mle_pdf_1_A30));
G_ub = Gv(find(cumsum(cond_mle_pdf_1_A30)>0.95,1,"first"));

display(['With 30% ascertainment and conditional index not ascertained, the Generation of detection: ' num2str(G_mle,'%d') ' (95 percentile = ' num2str(G_95,'%d') ')']);


% Generation: 50% Ascertainment 
Gv=[1:10];
G_mle = Gv(cond_mle_pdf_1_A50==max(cond_mle_pdf_1_A50));
G_95 = Gv(find(cumsum(cond_mle_pdf_1_A50)>0.95,1,"first"));

display(['With 50% ascertainment and conditional index not ascertained, the Generation of detection: ' num2str(G_mle,'%d') ' (95 percentile = ' num2str(G_95,'%d') ')']);

% Generation: 70% Ascertainment 
Gv=[1:10];
G_mle = Gv(cond_mle_pdf_1_A70==max(cond_mle_pdf_1_A70));
G_95 = Gv(find(cumsum(cond_mle_pdf_1_A70)>0.95,1,"first"));

display(['With 70% ascertainment and conditional index not ascertained, the Generation of detection: ' num2str(G_mle,'%d') ' (95 percentile = ' num2str(G_95,'%d') ')']);


% Generation: 90% Ascertainment 

Gv=[1:10];
G_mle = Gv(cond_mle_pdf_1_A90==max(cond_mle_pdf_1_A90));
G_95 = Gv(find(cumsum(cond_mle_pdf_1_A90)>0.95,1,"first"));

display(['With 90% ascertainment and conditional index not ascertained, the Generation of detection: ' num2str(G_mle,'%d') ' (95 percentile = ' num2str(G_95,'%d') ')']);
