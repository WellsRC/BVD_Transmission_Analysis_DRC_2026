clear;
clc;
load("BVD_Model_Fit.mat",'Parameters')
k=Parameters(:,1);
p=Parameters(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probability of extinction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prob_True_Extinction=p.^k;

Prob_False_Extinction = zeros(1001,length(k));
Prob_Apparent_Extinction = zeros(1001,length(k));
Conditional_Prob_Extinction = zeros(1001,length(k));
Odds_True_Extinction = zeros(1001,length(k));
uv=linspace(0,1,1001);
for ss=1:1001
    u=uv(ss);
    Prob_False_Extinction(ss,:)=p.^k.*(1-(1-u.*(1-p)).^k)./((1-u.*(1-p)).^k); 
    Prob_Apparent_Extinction(ss,:) = Prob_True_Extinction' + Prob_False_Extinction(ss,:);
    Conditional_Prob_Extinction(ss,:) = (1-u.*(1-p)).^k;
    Odds_True_Extinction(ss,:) = Prob_True_Extinction'.*(1-Prob_False_Extinction(ss,:))./((1-Prob_True_Extinction').*Prob_False_Extinction(ss,:));
end
close all;

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])


subplot('Position',[0.07 0.58 0.39 0.37])
lb = zeros(1,1001);
ub = zeros(1,1001);
mle = zeros(1,1001);

for ss=1:1001
    if(length(unique(Prob_False_Extinction(ss,:)))>1)
        bp = betafit(Prob_False_Extinction(ss,:));
        mle(ss) = (bp(1)-1)./(sum(bp)-2);
        lb(ss) = betainv(0.025,bp(1),bp(2));
        ub(ss) = betainv(0.975,bp(1),bp(2));
    else
        mle(ss) = unique(Prob_False_Extinction(ss,:));
        lb(ss) = unique(Prob_False_Extinction(ss,:));
        ub(ss) = unique(Prob_False_Extinction(ss,:));
    end
end

patch([uv flip(uv)], [lb flip(ub)],'k','LineStyle','none','FaceAlpha',0.3); hold on
plot(uv,mle,'k','LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:0.1:1],'Xminortick','on','Fontsize',16)

box off;
xlabel('Underascertainment rate','FontSize',14)
ylabel('Probability of false extinction','FontSize',14)
text(-0.18,1.02,'A','FontSize',28,'Units','normalized')

subplot('Position',[0.54 0.58 0.39 0.37])

lb = zeros(1,1001);
ub = zeros(1,1001);
mle = zeros(1,1001);

for ss=1:1000
    if(length(unique(Prob_Apparent_Extinction(ss,:)))>1)
        bp = betafit(Prob_Apparent_Extinction(ss,:));
        mle(ss) = (bp(1)-1)./(sum(bp)-2);
        lb(ss) = betainv(0.025,bp(1),bp(2));
        ub(ss) = betainv(0.975,bp(1),bp(2));
    else
        mle(ss) = unique(Prob_Apparent_Extinction(ss,:));
        lb(ss) = unique(Prob_Apparent_Extinction(ss,:));
        ub(ss) = unique(Prob_Apparent_Extinction(ss,:));
    end
end
ss=1001;
mle(ss) = 1;
        lb(ss) = 1;
        ub(ss) = 1;

patch([uv flip(uv)], [lb flip(ub)],'k','LineStyle','none','FaceAlpha',0.3); hold on
plot(uv,mle,'k','LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:0.1:1],'Xminortick','on','Fontsize',16)

box off;
xlabel('Underascertainment rate','FontSize',14)
ylabel('Probability of apparent extinction','FontSize',14)

text(-0.19,1.02,'B','FontSize',28,'Units','normalized')

subplot('Position',[0.07 0.09 0.39 0.37])

lb = zeros(1,1001);
ub = zeros(1,1001);
mle = zeros(1,1001);

for ss=1:1001
    if(length(unique(Conditional_Prob_Extinction(ss,:)))>1)
        bp = betafit(Conditional_Prob_Extinction(ss,:));
        mle(ss) = (bp(1)-1)./(sum(bp)-2);
        lb(ss) = betainv(0.025,bp(1),bp(2));
        ub(ss) = betainv(0.975,bp(1),bp(2));
    else
        mle(ss) = unique(Conditional_Prob_Extinction(ss,:));
        lb(ss) = unique(Conditional_Prob_Extinction(ss,:));
        ub(ss) = unique(Conditional_Prob_Extinction(ss,:));
    end
end

patch([uv flip(uv)], [lb flip(ub)],'k','LineStyle','none','FaceAlpha',0.3); hold on
plot(uv,mle,'k','LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:0.1:1],'Xminortick','on','Fontsize',16)

box off;
xlabel('Underascertainment rate','FontSize',14)
ylabel({'Probability of true extinction','given apparent extinction'},'FontSize',14)

text(-0.18,1.02,'C','FontSize',28,'Units','normalized')

subplot('Position',[0.54 0.09 0.39 0.37])

lb = zeros(1,1001);
ub = zeros(1,1001);
mle = zeros(1,1001);

for ss=1:1001
    if(length(unique(Odds_True_Extinction(ss,:)))>1)
        bp = gamfit(Odds_True_Extinction(ss,:));
        mle(ss) = (bp(1)-1).*bp(2);
        lb(ss) = gaminv(0.025,bp(1),bp(2));
        ub(ss) = gaminv(0.975,bp(1),bp(2));
    else
        mle(ss) = unique(Odds_True_Extinction(ss,:));
        lb(ss) = unique(Odds_True_Extinction(ss,:));
        ub(ss) = unique(Odds_True_Extinction(ss,:));
    end
end


patch([uv(2:end) flip(uv(2:end))], [lb(2:end) flip(ub(2:end))],'k','LineStyle','none','FaceAlpha',0.3); hold on
plot(uv(2:end),mle(2:end),'k','LineWidth',2);
plot(uv,ones(size(uv)),'-.','color',[0.5 0.5 0.5],'LineWidth',2);
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:0.1:1],'YTick',10.^[-1:4],'Xminortick','on','Fontsize',16)
yscale('log')
ylim([0.1 5.*10^4])
box off;
xlabel('Underascertainment rate','FontSize',14)
ylabel({'Odds of true extinction','relative to false extinction'},'FontSize',14)

text(-0.19,1.02,'D','FontSize',28,'Units','normalized')

theme('light')

print(gcf, 'Supplementary_Figure_2.png','-dpng','-r300');
