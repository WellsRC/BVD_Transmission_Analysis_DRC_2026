function BVD_Probability_of_False_Extinction(u_baseline,Scenario)
clc;
close all;


N_Samp=10^4;
[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p_base,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);

if(strcmp(Scenario,'Baseline'))
    p=p_base;
elseif(strcmp(Scenario,'Fatal'))
    p=p_fatal;
elseif(strcmp(Scenario,'Non-fatal'))
    p=p_non_fatal;
end
u=u_baseline.*ones(N_Samp,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comput the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prob_Extinction_Index=p.^k;

Prob_Asymptomatic_Transmission=p.^k.*(1-(1-u.*(1-p)).^k)./((1-u.*(1-p)).^k); 

A1=Prob_Extinction_Index./(Prob_Asymptomatic_Transmission+Prob_Extinction_Index);


N_Spillover=10;

p_temp=zeros(N_Samp,N_Spillover+1);

Events_Num=[0:N_Spillover];
for ii=1:N_Samp
    p_temp(ii,:)=binopdf(Events_Num,N_Spillover,1-A1(ii)); % Want to look at the numbr of potential asymptoamtic first generations
end

clc;

display(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
display([Scenario ': Under-ascertainment = ' num2str(u_baseline*100) '%']);
display(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
x=linspace(0,1,1001);

figure('units','normalized','outerposition',[0.05 0.05 0.9 1]);


subplot('Position',[0.05 0.575 0.44 0.41])

par_beta=betafit(Prob_Extinction_Index);

mle_p=(par_beta(1)-1)./(sum(par_beta)-2);
bnds_p=betainv([0.025 0.975],par_beta(1),par_beta(2));

display(['Probability of extinction at the index case: ' num2str(mle_p,'%4.3f') ' (95% UI: ' num2str(bnds_p(1),'%4.3f') char(8211) num2str(bnds_p(2),'%4.3f') ')'])

plot(x,betapdf(x,par_beta(1),par_beta(2)),'k','LineWidth',2);

set(gca,'LineWidth',2,'Fontsize',18,'tickdir','out')
ylabel('Probability density','FontSize',18)
xlabel('Probability of true extinction at index case','FontSize',18)
box off;
xlim([betainv([0.005 0.995],par_beta(1),par_beta(2))]);
text(-0.117,0.99,'A','FontSize',32,'VerticalAlignment','middle','Units','normalized');
subplot('Position',[0.55 0.575 0.44 0.41])

par_beta=betafit(Prob_Asymptomatic_Transmission);

mle_p=(par_beta(1)-1)./(sum(par_beta)-2);
bnds_p=betainv([0.025 0.975],par_beta(1),par_beta(2));

display(['Probability of false extinction: ' num2str(mle_p,'%4.3f') ' (95% UI: ' num2str(bnds_p(1),'%4.3f') char(8211) num2str(bnds_p(2),'%4.3f') ')'])


plot(x,betapdf(x,par_beta(1),par_beta(2)),'k','LineWidth',2);


set(gca,'LineWidth',2,'Fontsize',18,'tickdir','out');
ylabel('Probability density','FontSize',18)
xlabel('Probability of false extinction','FontSize',18)
box off;
xlim(betainv([0.0001 0.9999],par_beta(1),par_beta(2)));
text(-0.1255,0.99,'B','FontSize',32,'VerticalAlignment','middle','Units','normalized');

subplot('Position',[0.05 0.08 0.44 0.41])


par_beta=betafit(A1);

mle_p=(par_beta(1)-1)./(sum(par_beta)-2);
bnds_p=betainv([0.025 0.975],par_beta(1),par_beta(2));

display(['Likelihood of extinction relative to false-extinction: ' num2str(mle_p,'%4.3f') ' (95% UI: ' num2str(bnds_p(1),'%4.3f') char(8211) num2str(bnds_p(2),'%4.3f') ')'])


plot(x,betapdf(x,par_beta(1),par_beta(2)),'k','LineWidth',2);

set(gca,'LineWidth',2,'Fontsize',18,'tickdir','out');
ylabel('Probability density','FontSize',18)
xlabel('Likelihood of true extinction relative to false extinction','FontSize',18)
box off;
xlim(betainv([0.0001 0.9999],par_beta(1),par_beta(2)));
text(-0.117,0.99,'C','FontSize',32,'VerticalAlignment','middle','Units','normalized');

subplot('Position',[0.55 0.08 0.44 0.41])

beta_par=zeros(5,2); % Looking only at the first few number of chains that could produce unascertained first generations
for ii=1:5
    beta_par(ii,:)=betafit(p_temp(:,ii));
end

p_s=linspace(0.0001,0.9999,10^4);
L=zeros(10^4,1);
% Likelihood for first few number of chains that could produce unascertained first generations
for ii=1:10^4
    L(ii)=sum(log(betapdf(binopdf([0:4]',N_Spillover,p_s(ii)),beta_par(:,1),beta_par(:,2))));
end
q_mle=p_s(L==max(L));
w=exp(L)./sum(exp(L));
wc=cumsum(w);
Non_Extinct_Num=zeros(length(k),length(Events_Num));

for kk=1:length(k)
    r=rand(1);
    f_indx=find(r<=wc,1,'first');
    Non_Extinct_Num(kk,:)=binopdf(Events_Num,N_Spillover,p_s(f_indx));
end

lb_Non_Extinct_Num=binopdf(Events_Num,N_Spillover,q_mle)-prctile(Non_Extinct_Num,2.5,1);
ub_Non_Extinct_Num=prctile(Non_Extinct_Num,97.5,1)-binopdf(Events_Num,N_Spillover,q_mle);

bnds_p=[prctile(Non_Extinct_Num,[2.5 97.5],1)];

Non_Extinct_Num=binopdf(Events_Num,N_Spillover,q_mle);



mle_p=Non_Extinct_Num(1);




mle_c=Events_Num(Non_Extinct_Num==max(Non_Extinct_Num));
lb=Events_Num(find(cumsum(Non_Extinct_Num)<=0.025,1,'last'));
if(isempty(lb))
lb=0;
end
ub=Events_Num(find(cumsum(Non_Extinct_Num)>=0.975,1,'first'));
display(['FALSE extinction events among the ' num2str(N_Spillover) ' Cases: ' num2str(mle_c,'%1.0f') ' (95% CrI:' num2str(lb,'%1.0f') char(8211) num2str(ub,'%1.0f') ')' ]);

bar(Events_Num,Non_Extinct_Num,'k','LineStyle','none'); hold on;
er = errorbar(Events_Num,Non_Extinct_Num,lb_Non_Extinct_Num,ub_Non_Extinct_Num);    
er.Color = [0.7 0.7 0.7];                            
er.LineStyle = 'none';
er.LineWidth=2;

set(gca,'LineWidth',2,'Fontsize',18,'tickdir','out','YTick',[0:0.05:0.55]);
ylabel('Probability density','FontSize',18)
xlabel('Number of false exctintion events','FontSize',18)
box off;

xlim([-0.5 10.5])
ylim([0 0.55])

text(-0.1255,0.99,'D','FontSize',32,'VerticalAlignment','middle','Units','normalized');

theme('light')
print(gcf,['Supplementary_Figure_' Scenario '_Underascertainment=' num2str(u_baseline*100) '.png'],'-dpng','-r300');

end