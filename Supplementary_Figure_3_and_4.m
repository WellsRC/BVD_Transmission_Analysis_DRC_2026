clear;
clc;

% Ratio of the the likelihood of localized transmission relative to a
% spillover event

current_G=1;
x_past=[];

R0t=[0.05:0.05:2];
kt=10.^linspace(-1,1,101);
uv=0.01:0.01:0.99;

[k,u,R0]=meshgrid(kt,uv,R0t);
old_size=size(k);
k=k(:);
u=u(:);
R0=R0(:);
p=k./(k+R0);

Likelihood=zeros(length(k),1);
Relative_Risk=zeros(length(k),1);
Odds_Ratio=zeros(length(k),1);

for kk=1:length(k)
    zeta_G=zeros(50,1);
    for G=1:50
        if(G==1)
                zeta_G(G)=u(kk)-u(kk).*p(kk).^k(kk)./(1-u(kk).*(1-p(kk))).^k(kk);
        else
            temp_u = u(kk)-sum(zeta_G(1:G-1));
            zeta_G(G) = u(kk) -sum(zeta_G(1:G-1)) - u(kk).*p(kk).^k(kk)./(1-temp_u.*(1-p(kk))).^k(kk);
        end
    end

    L_HH=sum(zeta_G);
    L_S=1-u(kk);

    Likelihood(kk)=L_HH./(L_HH+L_S);
    Relative_Risk(kk)=L_HH./(L_S);
    Odds_Ratio(kk)=L_HH.*(1-L_S)./(L_S.*(1-L_HH));
end

cmap=[hex2rgb('#8e0152');
hex2rgb('#c51b7d');
hex2rgb('#de77ae');
hex2rgb('#f1b6da');
hex2rgb('#fde0ef');
hex2rgb('#f7f7f7');
hex2rgb('#e6f5d0');
hex2rgb('#b8e186');
hex2rgb('#7fbc41');
hex2rgb('#4d9221');
hex2rgb('#276419');];


cmap=flip(interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,1001)));

scaled_OR=Odds_Ratio;

scaled_OR=reshape(scaled_OR,old_size);

R0v=[0.5 1 1.5 2];
close all;
figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])

for rr=1:length(R0v)
    subplot('Position',[0.065+0.45.*rem((rr-1),2) 0.58-0.49.*floor((rr-1)/2) 0.39 0.37])
    temp_p=squeeze(scaled_OR(:,:,abs(R0t-R0v(rr))==min(abs(R0t-R0v(rr))))); % ran into numerical accuracy issues using R0t==R0v(rr)
    temp_p(temp_p>2)=2;
    contourf(uv,kt,temp_p',101,'LineStyle','none')
    yscale('log');
    clim([0 2]);
    colormap(cmap)
    hold on;


    [C,h] =contour(uv,kt,temp_p',[1],'k','LineWidth',2,'LineStyle',':');
    h.LevelList=1;
    clabel(C,h,'Fontsize',12,'LabelSpacing',1000,'Color','k');

    ylabel('Dispersion coefficient (k)','FontSize',14)
    xlabel('Underascertainment (\mu)','FontSize',14)
    set(gca,'Fontsize',12,'tickdir','out','XTick',[0.01 0.1:0.1:0.9 0.99],'LineWidth',2);
    yscale('log');
    box off;
    tp=title(['{\itR}_0 = ' num2str(R0v(rr),'%3.2f')]);
    tp.Position=tp.Position-[0 0.9 0];
    text(-0.125,1.02,char(65+rr-1),'FontSize',28,'Units','normalized');
end


subplot('Position',[0.93 0.09 0.015 0.86])

odds_HH=linspace(0,2,size(cmap,1));

for pp=1:length(odds_HH)-1
    patch([0 1 1 0],[odds_HH(pp) odds_HH(pp) odds_HH(pp+1) odds_HH(pp+1)],interp1(linspace(0,2,size(cmap,1)),cmap,odds_HH(pp)),'LineStyle','none');
end

odds_HH=[0:0.2:2];
for pp=1:length(odds_HH)-1
    text(1.2,odds_HH(pp),num2str(odds_HH(pp),'%3.2f'),'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
end

text(1.2,odds_HH(end),[num2str(odds_HH(end),'%3.2f') '+'],'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
text(3.5,1,'Odds ratio of localized transmission','Rotation',270,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
set(gca,'YTick',[],'XTick',[])
xlim([0 1])
ylim([0 2])


theme('light')
print(gcf,['Supplementary_Figure_3.png'],'-dpng','-r300');


pAv=[0.7 0.8 0.9 0.99];
close all;
figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])

for rr=1:length(pAv)
    subplot('Position',[0.065+0.45.*rem((rr-1),2) 0.58-0.49.*floor((rr-1)/2) 0.39 0.37])
    temp_p=squeeze(scaled_OR(abs(uv-pAv(rr))==min(abs(uv-pAv(rr))),:,:))'; % ran into numerical accuracy issues using R0t==R0v(rr)
    temp_p(temp_p>2)=2;
    contourf(R0t,kt,temp_p',101,'LineStyle','none')
    yscale('log');
    clim([0 2]);
    colormap(cmap)
    
    hold on;
    plot(ones(1,5),10.^linspace(-1,1,5),'Color',[0.75 0.75 0.75],'LineWidth',2,'LineStyle','-.')

    [C,h] =contour(R0t,kt,temp_p',[1],'k','LineWidth',2,'LineStyle',':');
    h.LevelList=1;
    clabel(C,h,'Fontsize',12,'LabelSpacing',1000,'Color','k');
    
    ylabel('Dispersion coefficient (k)','FontSize',14)
    xlabel('Basic reproduction number ({\itR}_0)','FontSize',14)
    set(gca,'Fontsize',12,'tickdir','out','XTick',[0.05 0.25:0.25:2],'LineWidth',2);
    yscale('log');
    box off;
    tp=title(['Underascertainment rate (\mu = ' num2str(pAv(rr),'%3.2f') ')']);
    tp.Position=tp.Position-[0 0.9 0];
    text(-0.125,1.02,char(65+rr-1),'FontSize',28,'Units','normalized');
end


subplot('Position',[0.93 0.09 0.015 0.86])

odds_HH=linspace(0,2,size(cmap,1));

for pp=1:length(odds_HH)-1
    patch([0 1 1 0],[odds_HH(pp) odds_HH(pp) odds_HH(pp+1) odds_HH(pp+1)],interp1(linspace(0,2,size(cmap,1)),cmap,odds_HH(pp)),'LineStyle','none');
end

odds_HH=[0:0.2:2];
for pp=1:length(odds_HH)-1
    text(1.2,odds_HH(pp),num2str(odds_HH(pp),'%3.2f'),'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
end

text(1.2,odds_HH(end),[num2str(odds_HH(end),'%3.2f') '+'],'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
text(3.5,1,'Odds ratio of localized transmission','Rotation',270,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
set(gca,'YTick',[],'XTick',[])
xlim([0 1])
ylim([0 2])
theme('light')
print(gcf,['Supplementary_Figure_4.png'],'-dpng','-r300');
