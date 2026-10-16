clear;
clc;
close all;


load('BVD_Model_Fit.mat')
R0 = Parameters(:,1).*(1-Parameters(:,2))./Parameters(:,2);
k= Parameters(:,2);
cfr = Parameters(:,3);
mu = Parameters(:,4);
I0 = Parameters(:,5);
G = Parameters(:,6);

figure('units','normalized','outerposition',[0.25 0.05 0.5 1])

subplot('Position',[0.1 0.715 0.4 0.265])
histogram(G,[-0.5:17.5],'FaceColor','k','FaceAlpha',1,'EdgeColor','w')
set(gca,'LineWidth',2,'tickdir','out','FontSize',12,'XTick',[0:17])
xlabel('Generation on May 13, 2026','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'A','FontSize',30,'Units','normalized')
box off;
xtickangle(0)
xlim([-0.5 17.5])

G_bin = [0:17];
N_temp = histcounts(G,[-0.5:17.5]);
mle = G_bin(N_temp==max(N_temp));
lb = prctile(G,2.5);
ub = prctile(G,97.5);

display(['Geneariton of May 13, 2026: ' num2str(mle,'%d') ' (95% CrI: ' num2str(lb,'%d') char(8211) num2str(ub,'%d') ')'])


subplot('Position',[0.58 0.715 0.4 0.265])
histogram(I0,[0.5:5.5],'FaceColor','k','FaceAlpha',1,'EdgeColor','w')
set(gca,'LineWidth',2,'tickdir','out','FontSize',12,'XTick',[1:5])
xlabel('Number of index cases','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'B','FontSize',30,'Units','normalized')
xlim([0.5 5.5])
box off;

I0_bin = [1:5];
N_temp = histcounts(I0,[0.5:5.5]);
mle = I0_bin(N_temp==max(N_temp));
lb = prctile(I0,2.5);
ub = prctile(I0,97.5);

display(['Number of index cases: ' num2str(mle,'%d') ' (95% CrI: ' num2str(lb,'%d') char(8211) num2str(ub,'%d') ')'])


subplot('Position',[0.1 0.38 0.4 0.265])
h=histfit(R0,51,"gamma");
set(gca,'LineWidth',2,'tickdir','out','FontSize',12,'XTick',[1.2:0.05:1.65])
xlabel('Basic reproduction number (\it{R}_0)','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'C','FontSize',30,'Units','normalized')
h(1).FaceColor=[0.5 0.5 0.5];
h(1).EdgeColor='w';
h(2).Color=hex2rgb('#C5001A');
h(2).LineStyle='-.';
box off;
xtickangle(0)
xlim([1.2 1.65])

par_dist = fitdist(R0,'Gamma');
mle = (par_dist.a-1).*par_dist.b;
lb = gaminv(0.025,par_dist.a,par_dist.b);
ub = gaminv(0.975,par_dist.a,par_dist.b);

display(['Basic reproduction number: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')'])


subplot('Position',[0.58 0.38 0.4 0.265])
h=histfit(k,51,"gamma");
set(gca,'LineWidth',2,'tickdir','out','FontSize',12,'XTick',[0.1:0.02:0.22])
xlabel('Dispersion coefficient (\it{k})','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'D','FontSize',30,'Units','normalized')
h(1).FaceColor=[0.5 0.5 0.5];
h(1).EdgeColor='w';
h(2).Color=hex2rgb('#C5001A');
h(2).LineStyle='-.';
box off;
xtickangle(0)
xlim([0.1 0.22])

par_dist = fitdist(k,'Gamma');
mle = (par_dist.a-1).*par_dist.b;
lb = gaminv(0.025,par_dist.a,par_dist.b);
ub = gaminv(0.975,par_dist.a,par_dist.b);

display(['Dispersion coefficient: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')'])

subplot('Position',[0.1 0.05 0.4 0.265])
h=histfit(1-mu,51,'kernel','Support',[0 1]);
set(gca,'LineWidth',2,'tickdir','out','FontSize',12,'XTick',[0.4:0.05:1])
xlabel('Ascertainment rate (1-\mu)','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'E','FontSize',30,'Units','normalized')
xlim([0.45 0.95])
h(1).FaceColor=[0.5 0.5 0.5];
h(1).EdgeColor='w';
h(2).Color=hex2rgb('#C5001A');
h(2).LineStyle='-.';
box off;
xtickangle(0)
% xlim([0 0.5])

par_dist = fitdist(1-mu,'kernel','Support',[0 1]);
mle = fmincon(@(z)-log(pdf(par_dist,z)),mean(1-mu),[],[],[],[],0,1);
lb = icdf(par_dist,0.025);
ub = icdf(par_dist,0.975);

display(['Ascertainment rate: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')'])

par_dist = fitdist(mu,'kernel','Support',[0 1]);
mle = fmincon(@(z)-log(pdf(par_dist,z)),mean(1-mu),[],[],[],[],0,1);
lb = icdf(par_dist,0.025);
ub = icdf(par_dist,0.975);

display(['Underascertainment rate: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')'])

subplot('Position',[0.58 0.05 0.4 0.265])
h=histfit(cfr,25,"beta");
set(gca,'LineWidth',2,'tickdir','out','FontSize',12)
xlabel('Case-fatality ratio (\chi)','FontSize',14)
ylabel('Frequency','FontSize',14)
text(-0.19,0.98,'F','FontSize',30,'Units','normalized')
h(1).FaceColor=[0.5 0.5 0.5];
h(1).EdgeColor='w';
h(2).Color=hex2rgb('#C5001A');
h(2).LineStyle='-.';
box off;
xtickangle(0)
xlim([0.3 0.6])

par_dist = fitdist(cfr,'beta');
mle = (par_dist.a-1)./(par_dist.a+par_dist.b-2);
lb = betainv(0.025,par_dist.a,par_dist.b);
ub = betainv(0.975,par_dist.a,par_dist.b);

display(['Case-fatality rate: ' num2str(mle,'%3.2f') ' (95% CrI: ' num2str(lb,'%3.2f') char(8211) num2str(ub,'%3.2f') ')'])


theme('light')

print(gcf,['Supplementary_Figure_1.png'],'-dpng','-r600');