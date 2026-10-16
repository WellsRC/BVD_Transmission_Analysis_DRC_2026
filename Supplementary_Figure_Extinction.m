clear;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comput the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0=[0.01:0.01:2];
k=10.^[-1:0.05:1];

Cmap=[hex2rgb('#ffffff');
    hex2rgb('#ffffe5');
hex2rgb('#f7fcb9');
hex2rgb('#d9f0a3');
hex2rgb('#addd8e');
hex2rgb('#78c679');
hex2rgb('#41ab5d');
hex2rgb('#238443');
hex2rgb('#006837');
hex2rgb('#004529');];

C_map=(interp1(linspace(0,1,size(Cmap,1)),Cmap,linspace(0,1,1001)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True extinction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])
Prob_Extinction=zeros(length(R0),length(k));
dE=zeros(3,length(R0),length(k));
Test=zeros(length(R0),length(k));
for gg=0:3
    for rr=1:length(R0)
        for kk=1:length(k)
            p=k(kk)./(R0(rr)+k(kk));
            if(gg==0)
                Prob_Extinction(rr,kk)=p.^k(kk);
            else
                dE(gg,rr,kk)=p.^k(kk)./((1-Prob_Extinction(rr,kk).*(1-p)).^k(kk))-Prob_Extinction(rr,kk);
                Prob_Extinction(rr,kk)=p.^k(kk)./((1-Prob_Extinction(rr,kk).*(1-p)).^k(kk));
            end
        end
    end
    subplot('Position',[0.065+0.45.*rem(gg,2) 0.58-0.49.*floor(gg/2) 0.39 0.37])
    contourf(R0,k,Prob_Extinction',101,'LineStyle','none');

    yscale('log');
    clim([0 1]);
    colormap(C_map)
    hold on;
    [C,h] =contour(R0,k,Prob_Extinction',[0:0.1:1],'k','LineWidth',2,'LineStyle',':');
     clabel(C,h,'Fontsize',12,'LabelSpacing',1000,'Color','k');
    plot(ones(1,5),10.^linspace(-1,1,5),'Color',[0.75 0.75 0.75],'LineWidth',2,'LineStyle','-.')
    ylabel('Dispersion coefficient (k)','FontSize',14)
    xlabel('Basic reproduction number ({\itR}_0)','FontSize',14)
    set(gca,'Fontsize',12,'tickdir','out','XTick',[0.01 0.25:0.25:2],'LineWidth',2);
    yscale('log');
    box off;
    switch gg
        case 0
            title(['Extinction at index case'])
        case 1
            title(['Extinction within one generation'])
        case 2
            title(['Extinction within two generations'])
        case 3
            title(['Extinction within three generations'])
    end
    text(-0.125,1.02,char(65+gg),'FontSize',28,'Units','normalized');
end

subplot('Position',[0.93 0.09 0.015 0.86])

p_ext=linspace(0,1,size(C_map,1));

for pp=1:length(p_ext)-1
    patch([0 1 1 0],[p_ext(pp) p_ext(pp) p_ext(pp+1) p_ext(pp+1)],interp1(linspace(0,1,size(Cmap,1)),Cmap,p_ext(pp)),'LineStyle','none');
end

p_ext=[0:0.1:1];
for pp=1:length(p_ext)
    text(1.2,p_ext(pp),num2str(p_ext(pp),'%3.2f'),'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
end
text(3.5,1.05/2,'Probability of extinction','Rotation',270,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
set(gca,'YTick',[],'XTick',[])
xlim([0 1])
ylim([0 1])

theme('light')

for gg=1:3
    for rr=1:2
        fprintf(['Absolute increase in probability of extinction for R0=' num2str(0.9+0.2.*(rr-1)) ' from generation ' num2str(gg-1) ' to generation ' num2str(gg) ':' num2str(squeeze(dE(gg,R0==0.9+0.2.*(rr-1),k==1)),'%4.3f') '\n']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% False extinction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cmap=[hex2rgb('#ffffff');
    hex2rgb('#ffffcc');
hex2rgb('#ffeda0');
hex2rgb('#fed976');
hex2rgb('#feb24c');
hex2rgb('#fd8d3c');
hex2rgb('#fc4e2a');
hex2rgb('#e31a1c');
hex2rgb('#bd0026');
hex2rgb('#800026');
hex2rgb('#000000');];
cmax_value=0.7;
C_map=(interp1(linspace(0,cmax_value,size(Cmap,1)),Cmap,linspace(0,cmax_value,1001)));

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])
Prob_False_Extinction=zeros(length(R0),length(k));
dE=zeros(3,length(R0),length(k));
for gg=0:3
    u=0.3+gg.*0.2;
    for rr=1:length(R0)
        for kk=1:length(k)
            p=k(kk)./(R0(rr)+k(kk));
            Prob_False_Extinction(rr,kk)=p.^k(kk).*(1-(1-u.*(1-p)).^k(kk))./((1-u.*(1-p)).^k(kk)); 
        end
    end
    subplot('Position',[0.065+0.45.*rem(gg,2) 0.58-0.49.*floor(gg/2) 0.39 0.37])
    contourf(R0,k,Prob_False_Extinction',101,'LineStyle','none')
    yscale('log');
    colormap(C_map)
    clim([0 cmax_value]);
    hold on;
    [C,h] =contour(R0,k,Prob_False_Extinction',[0:0.1:cmax_value],'k','LineWidth',2,'LineStyle',':');
     clabel(C,h,'Fontsize',12,'LabelSpacing',1000,'Color','k');
    plot(ones(1,5),10.^linspace(-1,1,5),'Color',[0.75 0.75 0.75],'LineWidth',2,'LineStyle','-.')
    ylabel('Dispersion coefficient (k)','FontSize',14)
    xlabel('Basic reproduction number ({\itR}_0)','FontSize',14)
    set(gca,'Fontsize',12,'tickdir','out','XTick',[0.01 0.25:0.25:2],'LineWidth',2);
    yscale('log');
    box off;
    title(['Underascertainment rate: ' num2str(u*100) '%' ])
    text(-0.125,1.02,char(65+gg),'FontSize',28,'Units','normalized');
end

subplot('Position',[0.93 0.09 0.015 0.86])

p_ext=linspace(0,cmax_value,size(C_map,1));

for pp=1:length(p_ext)-1
    patch([0 1 1 0],[p_ext(pp) p_ext(pp) p_ext(pp+1) p_ext(pp+1)],interp1(linspace(0,cmax_value,size(Cmap,1)),Cmap,p_ext(pp)),'LineStyle','none');
end

p_ext=[0:0.1:cmax_value];
for pp=1:length(p_ext)
    text(1.2,p_ext(pp),num2str(p_ext(pp),'%3.2f'),'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
end
text(3.5,cmax_value/2,'Probability of false extinction','Rotation',270,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
set(gca,'YTick',[],'XTick',[])
xlim([0 1])
ylim([0 cmax_value])

theme('light')

for gg=1:3
    for rr=1:2
        fprintf(['Absolute increase in probability of extinction for R0=' num2str(0.9+0.2.*(rr-1)) ' from generation ' num2str(gg-1) ' to generation ' num2str(gg) ':' num2str(squeeze(dE(gg,R0==0.9+0.2.*(rr-1),k==1)),'%4.3f') '\n']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True Extinction relative to False extinction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cmap=[hex2rgb('#ffffff');
    hex2rgb('#ffffe5');
hex2rgb('#f7fcb9');
hex2rgb('#d9f0a3');
hex2rgb('#addd8e');
hex2rgb('#78c679');
hex2rgb('#41ab5d');
hex2rgb('#238443');
hex2rgb('#006837');
hex2rgb('#004529');];

cmax_value=1;
C_map=(interp1(linspace(0,cmax_value,size(Cmap,1)),Cmap,linspace(0,cmax_value,1001)));

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])
Prob_False_Extinction=zeros(length(R0),length(k));
Prob_True_Extinction=zeros(length(R0),length(k));
for gg=0:3
    u=0.3+gg.*0.2;
    for rr=1:length(R0)
        for kk=1:length(k)
            p=k(kk)./(R0(rr)+k(kk));
            Prob_False_Extinction(rr,kk)=p.^k(kk).*(1-(1-u.*(1-p)).^k(kk))./((1-u.*(1-p)).^k(kk)); 
            Prob_True_Extinction(rr,kk)=p.^k(kk);
        end
    end
    R=Prob_True_Extinction./(Prob_False_Extinction+Prob_True_Extinction);
    subplot('Position',[0.065+0.45.*rem(gg,2) 0.58-0.49.*floor(gg/2) 0.39 0.37])
    contourf(R0,k,R',101,'LineStyle','none')
    yscale('log');
    colormap(C_map)
    clim([0 cmax_value]);
    hold on;
    [C,h] =contour(R0,k,R',[0:0.1:cmax_value],'k','LineWidth',2,'LineStyle',':');
     clabel(C,h,'Fontsize',12,'LabelSpacing',1000,'Color','k');
    plot(ones(1,5),10.^linspace(-1,1,5),'Color',[0.75 0.75 0.75],'LineWidth',2,'LineStyle','-.')
    ylabel('Dispersion coefficient (k)','FontSize',14)
    xlabel('Basic reproduction number ({\itR}_0)','FontSize',14)
    set(gca,'Fontsize',12,'tickdir','out','XTick',[0.01 0.25:0.25:2],'LineWidth',2);
    yscale('log');
    box off;
    title(['Underascertainment rate: ' num2str(u*100) '%' ])
    text(-0.125,1.02,char(65+gg),'FontSize',28,'Units','normalized');
end

subplot('Position',[0.93 0.09 0.015 0.86])

p_ext=linspace(0,cmax_value,size(C_map,1));

for pp=1:length(p_ext)-1
    patch([0 1 1 0],[p_ext(pp) p_ext(pp) p_ext(pp+1) p_ext(pp+1)],interp1(linspace(0,cmax_value,size(Cmap,1)),Cmap,p_ext(pp)),'LineStyle','none');
end

p_ext=[0:0.1:cmax_value];
for pp=1:length(p_ext)
    text(1.2,p_ext(pp),num2str(p_ext(pp),'%3.2f'),'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle')
end
text(3.5,cmax_value/2,'Probability of true extinction relative to false extinction','Rotation',270,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
set(gca,'YTick',[],'XTick',[])
xlim([0 1])
ylim([0 cmax_value])

theme('light')



load("BVD_Model_Fit.mat",'Parameters')

g_k =gamfit(Parameters(:,1));
k_mle_base =(g_k(1)-1).*g_k(2);
k = Parameters(:,1);
p=Parameters(:,2);

Prob_False_Extinction=zeros(10^3,100);
Prob_True_Extinction=zeros(10^3,100);
for rr=1:10^3
    for uu=1:100
        u=(uu)./100;
        Prob_False_Extinction(rr,uu)=p(rr).^k(rr).*(1-(1-u.*(1-p(rr))).^k(rr))./((1-u.*(1-p(rr))).^k(rr)); 
        Prob_True_Extinction(rr,uu)=p(rr).^k(rr);
    end
end
R=Prob_True_Extinction.*(1-Prob_False_Extinction)./(Prob_False_Extinction.*(1-Prob_True_Extinction));

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])
subplot('Position',[0.07 0.10 0.9 0.88])
uv=linspace(0.01,1,100);

plot(uv,ones(size(uv)),':','LineWidth',1.5,'Color',[0.5 0.5 0.5]); hold on;

patch([uv flip(uv)],[prctile(R,2.5,1) flip(prctile(R,97.5,1))],'k','FaceAlpha',0.3,'linestyle','none'); hold on
patch([uv flip(uv)],[prctile(R,25,1) flip(prctile(R,75,1))],'k','FaceAlpha',0.3,'linestyle','none'); hold on
plot(uv,median(R),'k-','LineWidth',2)


ylabel('Odds of true extinction relative to false extinction','FontSize',14)
xlabel('Underascertainment rate (\mu)','FontSize',14)
set(gca,'Fontsize',12,'tickdir','out','XTick',[0.01 0.1:0.1:1],'LineWidth',2);
xlim([0.01 1])
ylim([10^(-1) 2.*10^3])
yscale('log');
box off;
theme('light')
