clear;
clc;
close all;

load('BVD_Model_Fit.mat')
Model_Cases=cumsum(Model_Cases,2);
Model_Deaths=cumsum(Model_Deaths,2);
Generation = Parameters(:,6);
figure('units','normalized','outerposition',[0.25 0.05 0.5 0.9])
subplot('Position',[0.1 0.56 0.87 0.42])
for tt=1:length(Date)    
    lb = prctile(Model_Cases(:,tt),2.5);
    ub = prctile(Model_Cases(:,tt),97.5);
   p1= patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#002C54'),'FaceAlpha',0.3,'LineStyle','none'); hold on

    lb = prctile(Model_Cases(:,tt),25);
    ub = prctile(Model_Cases(:,tt),75);
    patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#002C54'),'FaceAlpha',0.3,'LineStyle','none'); hold on


    lb = prctile(Model_Deaths(:,tt),2.5);
    ub = prctile(Model_Deaths(:,tt),97.5);
   p2 = patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#1E0000'),'FaceAlpha',0.3,'LineStyle','none'); hold on
    
    lb = prctile(Model_Deaths(:,tt),25);
    ub = prctile(Model_Deaths(:,tt),75);
    patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#1E0000'),'FaceAlpha',0.3,'LineStyle','none'); hold on

    
end

   scatter(1:length(Date),cumsum(Cases),30,hex2rgb('#16253D'),'filled');
   scatter(1:length(Date),cumsum(Deaths),30,hex2rgb('#500805'),'filled','Marker','diamond');
set(gca,'LineWidth',2,'TickDir','out','XTick',[1:length(Date)],'XTickLabel',{datestr(Date)},'FontSize',12,'YMinorTick','on')
ylabel('Cumulative ascertained incidence','FontSize',14)
xlabel('Date','FontSize',14)
legend([p1 p2],{'Cases','Deaths'},'FontSize',14,'Location','northwest')
legend boxoff
xlim([0.5 length(Date)+0.5])
subplot('Position',[0.1 0.07 0.87 0.42])

for tt=1:length(Date)
    lb = prctile(Model_Cases_Overall(:,tt)-Model_Cases(:,tt),2.5);
    ub = prctile(Model_Cases_Overall(:,tt)-Model_Cases(:,tt),97.5);
    patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#002C54'),'FaceAlpha',0.3,'LineStyle','none'); hold on

    lb = prctile(Model_Cases_Overall(:,tt)-Model_Cases(:,tt),25);
    ub = prctile(Model_Cases_Overall(:,tt)-Model_Cases(:,tt),75);
    patch(tt+[-0.4 -0.4 0.4 0.4],[lb ub ub lb],hex2rgb('#002C54'),'FaceAlpha',0.3,'LineStyle','none'); hold on
end
ax = gca; 
ax.YAxis.Exponent = 0;
set(gca,'LineWidth',2,'TickDir','out','XTick',[1:length(Date)],'XTickLabel',{datestr(Date)},'FontSize',12,'YMinorTick','on')
ylabel('Cumulative unascertained incidence','FontSize',14)
xlabel('Date','FontSize',14)
xlim([0.5 length(Date)+0.5])
subplot('Position',[0.1 0.07 0.87 0.42])
theme('light')

print(gcf,['Figure_1.png'],'-dpng','-r600');


gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,15.36/(z-1))-[13.29 17.42],100,1,5000);
a_gp_serial=gp;
b_gp_serial=15.36/(gp-1);
gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,7.68/(z-1))-[7.04 8.32],100,1,5000);
a_gp_incubation=gp;
b_gp_incubation=7.68/(gp-1);
gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,10.77/(z-1))-[9.74 11.80],100,1,5000);
a_gp_death=gp;
b_gp_death=10.77/(gp-1);

test_time=zeros(length(Generation),1);

for ii=1:length(Generation)

    d = gamrnd(a_gp_death,b_gp_death);
    serial_samp = gamrnd(a_gp_serial,b_gp_serial,1,Generation(ii));
    inc_samp = gamrnd(a_gp_incubation,b_gp_incubation,1);
    temp_samp =  inc_samp + sum(serial_samp);
    temp_samp_2 = inc_samp + sum(serial_samp(1:end-1));
    test_time(ii) =max(d + temp_samp_2,temp_samp);
end


[N,edges] = histcounts(test_time,-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated: ' datestr(mle) ' (95% CrI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);