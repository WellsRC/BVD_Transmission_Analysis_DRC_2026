clear;
clc;
if isempty(gcp('nocreate'))
    parpool(6);
end

load("BVD_Model_Fit.mat",'Parameters')
NS=size(Parameters,1);

R0=Parameters(:,1).*(1-Parameters(:,2))./Parameters(:,2);

k=Parameters(:,1);
p=Parameters(:,2);
u=Parameters(:,4);

delta_I=1/10;

gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,15.36/(z-1))-[13.29 17.42],100,1,5000);
a_gp_serial=gp;
b_gp_serial=15.36/(gp-1);
g=15.36; % Serial interval as a proxy for generation time

t_ref=linspace(0,400,401);
Prev = zeros(NS,length(t_ref));
I0=10^(-4);
S0=1-I0;
opts =odeset('RelTol',10^(-12),'AbsTol',10^(-12),'NonNegative',1:2);
for uu=1:NS
    [T,Y]=ode45(@(t,x) SIR_ODE(t,x,R0(uu),delta_I), [0 400], [S0 I0],opts);
    Prev(uu,:)=exp(pchip(T,log(Y(:,2)),t_ref));
end

b = betafit(u);
uv = (b(1)-1)./(sum(b)-2);

PDF_prev = zeros(1,length(t_ref));
PDF_prev_lb = zeros(1,length(t_ref));
PDF_prev_ub = zeros(1,length(t_ref));

Past_prev = zeros(NS,length(t_ref));

f_discrete=[7:7:28];

Cumulative_PDF_prev = zeros(4,length(t_ref));
Cumulative_PDF_prev_lb = zeros(4,length(t_ref));
Cumulative_PDF_prev_ub = zeros(4,length(t_ref));

Number_Unknown_Exposures=1;

parfor tt=1:length(t_ref)
    % Ratio of the the likelihood of human-to-human transmission relative to a
    % spillover event
    R=zeros(NS,1);

    current_G=1;
    x_past=[];
    for ss=1:NS
        zeta_G=zeros(10,1);
        past_seeding_rate=zeros(10,1);
        for G=1:10
            if(G==1)
                if(t_ref(tt)-g*G<0)
                    seed_rate=0;
                else    
                    seed_rate = pchip(t_ref,Prev(ss,:),t_ref(tt)-g.*G);
                end
                zeta_G(G)=seed_rate.*(u(ss)-u(ss).*p(ss).^k(ss)./(1-u(ss).*(1-p(ss))).^k(ss));
            else
                past_seeding_rate_temp = past_seeding_rate(1:(G-1));
                zeta_ratio=zeta_G(1:G-1)./past_seeding_rate_temp(:);
                zeta_ratio(past_seeding_rate_temp(:)==0)=0;
                temp_u = u(ss)-sum(zeta_ratio);

                if(t_ref(tt)-g*G<0)
                    seed_rate=0;
                else    
                    seed_rate = pchip(t_ref,Prev(ss,:),t_ref(tt)-g.*G);
                end

                zeta_G(G) = seed_rate.*(temp_u - u(ss).*p(ss).^k(ss)./(1-temp_u.*(1-p(ss))).^k(ss));
            end
            past_seeding_rate(G)=seed_rate;
        end

        L_HH=sum(zeta_G);
        L_S=(1-u(ss)).*pchip(t_ref,Prev(ss,:),t_ref(tt));

        R(ss)=L_HH./(L_HH+L_S);
    end

    temp_p=1-(1-R).^Number_Unknown_Exposures;


    Past_prev(:,tt)=temp_p;



    temp_p(temp_p==1)=1-10^(-10);
    if(length(unique(temp_p))>1)
        beta_HH_alt=betafit(temp_p);
        PDF_prev(tt)=(beta_HH_alt(1)-1)./(sum(beta_HH_alt)-2);
        PDF_prev_lb(tt)=betainv(0.025,beta_HH_alt(1),beta_HH_alt(2));
        PDF_prev_ub(tt)=betainv(0.975,beta_HH_alt(1),beta_HH_alt(2));
    else
        PDF_prev(tt)=unique(temp_p);
        PDF_prev_lb(tt)=unique(temp_p);
        PDF_prev_ub(tt)=unique(temp_p);
    end
end

for ww=1:length(f_discrete)
    X=zeros(1,length(t_ref));
    for ii=1:401
        X= 1-prod((1-Past_prev(:,1:ii)./f_discrete(ww)),2);
        if(length(unique(X))>1)
                z = log(X./(1-X));
                [m_temp,s_temp] = normfit(z);
                Cumulative_PDF_prev(ww,ii) =1./(1+exp(-m_temp));
                l_temp=norminv(0.025,m_temp,s_temp);
                u_temp=norminv(0.975,m_temp,s_temp);
                Cumulative_PDF_prev_lb(ww,ii)=1./(1+exp(-l_temp));
                Cumulative_PDF_prev_ub(ww,ii)=1./(1+exp(-u_temp));
        else
            Cumulative_PDF_prev(ww,ii) =unique(X);
            Cumulative_PDF_prev_lb(ww,ii)=unique(X);
            Cumulative_PDF_prev_ub(ww,ii)=unique(X);
        end
    end
end

display(['Probability local transmission vs seeding event (EARLY): ' num2str(PDF_prev(t_ref==16),'%4.3f') ' (95% CrI: ' num2str(PDF_prev_lb(t_ref==16),'%4.3f') char(8211) num2str(PDF_prev_ub(t_ref==16),'%4.3f') ')'])
display(['Probability local transmission vs seeding event (LATE): ' num2str(PDF_prev(t_ref==400),'%4.3f') ' (95% CrI: ' num2str(PDF_prev_lb(t_ref==400),'%4.3f') char(8211) num2str(PDF_prev_ub(t_ref==400),'%4.3f') ')'])

for ww=1:length(f_discrete)
    display(['(Frequency ' num2str(f_discrete(ww),'%d') ') Cumualtive Probability local transmission vs seeding event at t=195: ' num2str(Cumulative_PDF_prev(ww,t_ref==195),'%4.3f') ' (95% CrI: ' num2str(Cumulative_PDF_prev_lb(ww,t_ref==195),'%4.3f') char(8211) num2str(Cumulative_PDF_prev_ub(ww,t_ref==195),'%4.3f') ')'])
end


close all;
figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])

subplot('Position',[0.065 0.58 0.38 0.37])
patch([t_ref flip(t_ref)],[prctile(Prev,[25]) flip(prctile(Prev,[75]))],'k','LineStyle','none','FaceAlpha',0.3); hold on
l2=patch([t_ref flip(t_ref)],[prctile(Prev,[2.5]) flip(prctile(Prev,[97.5]))],'k','LineStyle','none','FaceAlpha',0.3); hold on;

l1=patch([-10 -10 -9 -9],[-10 -9 -9 -10],'k','LineStyle','none','FaceAlpha',0.6);
l3=plot(t_ref,mean(Prev),'k','LineWidth',2); hold on;
l4=plot(t_ref,median(Prev),'k-.','LineWidth',2); hold on
set(gca,'TickDir','out','LineWidth',2,'XTick',[0:50:400],'YTick',[0:0.01:0.1],'Fontsize',14,'Xminortick','on')
xlim([0 400])
ylim([0 0.08])
ylabel({'Prevalence at epicenter'},'Fontsize',16)
xlabel('Day','Fontsize',16)
box off
legend([l3 l4 l2 l1],{'Mean','Median','95% CrI','IQR'})
legend boxoff;
xtickangle(0)

text(-0.16,0.99,'A','FontSize',28,'VerticalAlignment','middle','Units','normalized');

subplot('Position',[0.525 0.58 0.38 0.37])
patch([t_ref flip(t_ref)],[PDF_prev_lb flip(PDF_prev_ub)],'k','LineStyle','none','FaceAlpha',0.3); hold on
plot(t_ref,PDF_prev,'color','k','LineWidth',2); hold on
set(gca,'TickDir','out','LineWidth',2,'XTick',[0:50:400],'YTick',[0:0.05:0.25],'Fontsize',14,'Xminortick','on')
xlim([0 400])
ylim([0 0.25])
ylabel({'Probability',' of local transmission'},'Fontsize',16)
xlabel('Day','Fontsize',16)
box off
text(-0.2,0.99,'B','FontSize',28,'VerticalAlignment','middle','Units','normalized');
xtickangle(0)


subplot('Position',[0.525 0.09 0.38 0.37])

f =linspace(7,35,29);
Y_PDF_Normalized = zeros(length(f),length(t_ref));
T_MAX = zeros(length(f),1);
for ww=1:length(f)
    X=zeros(NS,length(t_ref));
    for ii=1:401
        X(:,ii)= 1-prod((1-Past_prev(:,1:ii)./f(ww)),2);
    end
    z = log(X./(1-X));
    [m_temp,s_temp] = normfit(z);
    X_mle =1./(1+exp(-m_temp));
    Y_PDF_Normalized(ww,:) = [X_mle(1) diff(X_mle)];
     Y_PDF_Normalized(ww,:)= Y_PDF_Normalized(ww,:)./max( Y_PDF_Normalized(ww,:));
     T_MAX(ww)= t_ref(Y_PDF_Normalized(ww,:)==max(Y_PDF_Normalized(ww,:)));
end
test = max(Y_PDF_Normalized,[],2);
pcolor(t_ref,f,Y_PDF_Normalized,'LineStyle','none')

c=colormap("pink");
colormap(flip(c));
c= colorbar;
c.Position=[0.912846325779309,0.08994708994709,0.015873015873016,0.37037037037037];

c.Label.String ={"Normalized probability","of local transmission"};
c.Label.Rotation=270;
c.Label.FontSize=16;

set(gca,'TickDir','out','LineWidth',2,'XTick',[0:50:400],'YTick',[7:7:35],'Yminortick','on','Fontsize',14,'Xminortick','on')
xlim([0 400])
box off;
xtickangle(0)

ylabel({'Frequency of cases','with unknown exposure'},'Fontsize',16)
xlabel('Day','Fontsize',16)

text(-0.2,0.99,'D','FontSize',28,'VerticalAlignment','middle','Units','normalized');
subplot('Position',[0.065 0.09 0.38 0.37])

plot(t_ref,Cumulative_PDF_prev(f_discrete==7,:),'Color',hex2rgb('#444C5C'),'LineWidth',2); hold on
plot(t_ref,Cumulative_PDF_prev(f_discrete==14,:),'Color',hex2rgb('#CE5A57'),'LineWidth',2)
plot(t_ref,Cumulative_PDF_prev(f_discrete==21,:),'Color',hex2rgb('#78A5A3'),'LineWidth',2)
plot(t_ref,Cumulative_PDF_prev(f_discrete==28,:),'Color',hex2rgb('#E1B16A'),'LineWidth',2)

set(gca,'TickDir','out','LineWidth',2,'XTick',[0:50:400],'YTick',[0:0.1:1],'Fontsize',14,'Xminortick','on')
ylabel({'Cumulative probability','of local transmission'},'Fontsize',14)
xlabel('Day','Fontsize',16)
box off;
ylim([0 1])
legend({'7-days','14-days','21-days','28-days'},'Location','northwest')
legend boxoff;
text(-0.16,0.99,'C','FontSize',28,'VerticalAlignment','middle','Units','normalized');
theme('light')
% 
% 
print(gcf,['Figure_2.png'],'-dpng','-r600');