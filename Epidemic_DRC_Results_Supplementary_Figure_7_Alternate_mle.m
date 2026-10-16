clear;

load('BVD_Model_Fit.mat','Parameters');

k=Parameters(:,1);
p=Parameters(:,2);
cfr = Parameters(:,3);
mu =Parameters(:,4);
g_k =gamfit(Parameters(:,1));
k_mle=(g_k(1)-1).*g_k(2);
beta_p = betafit(Parameters(:,2));
p_mle = (beta_p(1)-1)./(sum(beta_p)-2);
beta_p = betafit(Parameters(:,3));
cfr_mle = (beta_p(1)-1)./(sum(beta_p)-2);
beta_p = betafit(Parameters(:,4));
mu_mle = (beta_p(1)-1)./(sum(beta_p)-2);
N_Samp=10^4;


% Construct a proxy incidence curve from the historical information 
rng(20260720)
for IndxC=1:5
    for ss=1:IndxC
        if(ss==1)
            [Incidence] = Simulate_Branching_Epidemic(k_mle.*ones(N_Samp,1),p_mle.*ones(N_Samp,1),20);
        else            
            Incidence = Incidence+Simulate_Branching_Epidemic(k_mle.*ones(N_Samp,1),p_mle.*ones(N_Samp,1),20);
        end
    end
    if(IndxC==1)
        Incidence_Prev_Approx = Incidence;
    end
end

load('Epidemic_Analysis_Likelihood_Alternate_mle','Lt_Lagged','Lt_Generation');

Indx=[1:5];
Generation =[0:15];
[Generation,Indx]=meshgrid(Generation,Indx);
Indx=Indx(:);
Generation=Generation(:);
% underascertain_rate=underascertain_rate(:);
% prob_underascertain=prob_underascertain(:);

% samp_underascertain_rate=zeros(10^6,3);
% samp_prob_underascertain=zeros(10^6,3);
samp_Indx=zeros(10^6,3);
samp_Generation=zeros(10^6,3);
w=Lt_Lagged(:)./sum(Lt_Lagged(:));
wc=cumsum(w);
for ii=1:10^6
    f_indx = find(rand(1)<=wc,1,"first");
    % samp_underascertain_rate(ii,1)=underascertain_rate(f_indx);
    % samp_prob_underascertain(ii,1)=prob_underascertain(f_indx);
    samp_Indx(ii,1)=Indx(f_indx);
    samp_Generation(ii,1)=Generation(f_indx);
end

w=Lt_Generation(:)./sum(Lt_Generation(:));
wc=cumsum(w);
for ii=1:10^6
    f_indx = find(rand(1)<=wc,1,"first");
    % samp_underascertain_rate(ii,2)=underascertain_rate(f_indx);
    % samp_prob_underascertain(ii,2)=prob_underascertain(f_indx);
    samp_Indx(ii,2)=Indx(f_indx);
    samp_Generation(ii,2)=Generation(f_indx);
end


w=Lt_Generation(:).*Lt_Lagged(:)./sum(Lt_Generation(:).*Lt_Lagged(:));
wc=cumsum(w);
for ii=1:10^6
    f_indx = find(rand(1)<=wc,1,"first");
    % samp_underascertain_rate(ii,3)=underascertain_rate(f_indx);
    % samp_prob_underascertain(ii,3)=prob_underascertain(f_indx);
    samp_Indx(ii,3)=Indx(f_indx);
    samp_Generation(ii,3)=Generation(f_indx);
end

[N,edges] = histcounts(samp_Generation(:,1),-0.5:15.5);
G_f = edges(1:end-1)+0.5;
G_ref = G_f(N==max(N));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRC Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=zeros(N_Samp,1); % For the probability of local transmission comapred to seeding event
Prev_Approx=Incidence_Prev_Approx; % Use incidence as proxy for prevalence
c_temp =cumsum(Incidence_Prev_Approx,2); 

indx_shift_generation= 1; % The shift  G_ref to add the additional generation/serial interval
G_ref_DRC = G_ref+indx_shift_generation; % May 24 would be the next generation from the referenced
f_indx = c_temp(:,G_ref+1)>=246; % Find the simulations where there are at least 246 cases at the reference generation (add 1 because the index Gref could be zero)
Prev_Approx = Prev_Approx(f_indx,:); % Set the prevalence to all of those scenarioes
u=0.57; % The underascertainment rate
for ss=1:N_Samp
    zeta_G=zeros(10,1); % Probability of local transmission based on generation
    past_seeding_rate=zeros(10,1); % past seeding rates
    for G=1:10 % Numebr of generations back
        if(G==1)
            if(G_ref_DRC-G<0)
                seed_rate=0;
            elseif(G_ref_DRC==G)
                seed_rate=mean(Prev_Approx(:,G_ref_DRC-G+1)); % Initial condition (i.e. generation zero)
            else
                seed_rate = Prev_Approx(randi(nnz(f_indx)),G_ref_DRC-G+1);
            end
            zeta_G(G)=seed_rate.*(u-u.*p_mle.^k_mle./(1-u.*(1-p_mle)).^k_mle);
        else
            past_seeding_rate_temp = past_seeding_rate(1:(G-1));
            zeta_ratio=zeta_G(1:G-1)./past_seeding_rate_temp(:);
            zeta_ratio(past_seeding_rate_temp(:)==0)=0;
            temp_u = u-sum(zeta_ratio);

            if(G_ref_DRC-G<0)
                seed_rate=0;
            elseif(G_ref_DRC==G)
                seed_rate=mean(Prev_Approx(:,G_ref_DRC-G+1)); % Initial condition
            else
                seed_rate = Prev_Approx(randi(nnz(f_indx)),G_ref_DRC-G+1);
            end

            zeta_G(G) = seed_rate.*(temp_u - u.*p_mle.^k_mle./(1-temp_u.*(1-p_mle)).^k_mle);
        end
        past_seeding_rate(G)=seed_rate;
    end

    L_HH=sum(zeta_G);
    L_S=(1-u).*Prev_Approx(randi(nnz(f_indx)),G_ref_DRC+1);

    R(ss)=L_HH./(L_HH+L_S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uganda Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://academic.oup.com/cid/article-abstract/62/10/1277/2462818?login=false
u_mle =0.271;
bp_asym = lsqnonlin(@(z) betainv([0.025 0.975],z,(z-u_mle.*z+2.*u_mle-1)./u_mle)-[0.145 0.396],3,1,100);
b=(bp_asym-u_mle.*bp_asym+2.*u_mle-1)./u_mle;
u=betarnd(bp_asym,b,N_Samp,1);
indx_shift_generation= 1; % The shift  G_ref to add the additional generation/serial interval
G_ref_U = G_ref+indx_shift_generation; % May 24 would be the next generation from the reference generation for May 13

x_past=[];
R_Uganda = zeros(N_Samp,1);
Prev_Approx=Incidence_Prev_Approx;
c_temp =cumsum(Incidence_Prev_Approx,2);
Prev_Approx = Prev_Approx(f_indx,:);
Cases_Uganda_Reported=7;  % https://reliefweb.int/report/democratic-republic-congo/unicef-ebola-virus-disease-outbreak-flash-update-no-1-democratic-republic-congo-uganda-26-may-2026
% Cases_Uganda_Reported=2; %% https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON603  (May 21)
Cases_Uganda_Local= zeros(length(R_Uganda),10^3);
for ss=1:N_Samp
    zeta_G=zeros(10,1);
    past_seeding_rate=zeros(10,1);
    for G=1:10 % Numebr of generations back
        if(G==1)
            if(G_ref_U-G<0)
                seed_rate=0;
            elseif(G_ref_U==G)
                seed_rate=mean(Prev_Approx(:,G_ref_U-G+1)); % Initial condition
            else
                seed_rate = Prev_Approx(randi(nnz(f_indx)),G_ref_U-G+1);
            end
            zeta_G(G)=seed_rate.*(u(ss)-u(ss).*p_mle.^k_mle./(1-u(ss).*(1-p_mle)).^k_mle);
        else
            past_seeding_rate_temp = past_seeding_rate(1:(G-1));
            zeta_ratio=zeta_G(1:G-1)./past_seeding_rate_temp(:);
            zeta_ratio(past_seeding_rate_temp(:)==0)=0;
            temp_u = u(ss)-sum(zeta_ratio);

            if(G_ref_U-G<0)
                seed_rate=0;
            elseif(G_ref_U==G)
                seed_rate=mean(Prev_Approx(:,G_ref_U-G+1)); % Initial condition
            else
                seed_rate = Prev_Approx(randi(nnz(f_indx)),G_ref_U-G+1);
            end

            zeta_G(G) = seed_rate.*(temp_u - u(ss).*p_mle.^k_mle./(1-temp_u.*(1-p_mle)).^k_mle);
        end
        past_seeding_rate(G)=seed_rate;
    end

    L_HH=sum(zeta_G);
    L_S=(1-u(ss)).*Prev_Approx(randi(nnz(f_indx)),G_ref_U+1);

    R_Uganda(ss)=L_HH./(L_HH+L_S);
       Cases_Uganda_Local(ss,:)=binornd(Cases_Uganda_Reported,R_Uganda(ss),1,10^3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time distributions for determining the origin time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,15.36/(z-1))-[13.29 17.42],100,1,5000);
a_gp_serial=gp;
b_gp_serial=15.36/(gp-1);
gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,7.68/(z-1))-[7.04 8.32],100,1,5000);
a_gp_incubation=gp;
b_gp_incubation=7.68/(gp-1);
gp = lsqnonlin(@(z) gaminv([0.025 0.975],z,10.77/(z-1))-[9.74 11.80],100,1,5000);
a_gp_death=gp;
b_gp_death=10.77/(gp-1);

test_time=zeros(length(samp_Generation(:,1)),3);
test_time_alt=zeros(length(samp_Generation(:,1)),3);
for ii=1:length(samp_Generation(:,1))
    d = gamrnd(a_gp_death,b_gp_death);
    for mm=1:3
        serial_samp = gamrnd(a_gp_serial,b_gp_serial,1,samp_Generation(ii,mm));
        inc_samp = gamrnd(a_gp_incubation,b_gp_incubation,1);
        temp_samp =  inc_samp + sum(serial_samp);
        temp_samp_2 = inc_samp + sum(serial_samp(1:end-1));
        test_time(ii,mm) = max(d + temp_samp_2,temp_samp);
    end
end

clc;
clc;
display('=================================================================================')
display('Lagged')
display('=================================================================================')
[N,edges] = histcounts(samp_Indx(:,1),0.5:5.5);
mle = find(N==max(N));
test_bnds = cumsum(N)./sum(N);
bnds = [find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")];
display(['Number of index cases: ' num2str(mle,'%d') ' (95% CrI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);

[N,edges] = histcounts(samp_Generation(:,1),-0.5:15.5);
G_f = edges(1:end-1)+0.5;
mle = G_f(find(N==max(N)));
mle_Gen = mle;
test_bnds = cumsum(N)./sum(N);
bnds = G_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]);
display(['Generation of epidemic: ' num2str(mle,'%d') ' (95% CrI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);

bp = betafit(R(:));
mle = (bp(1)-1)./(sum(bp)-2);
bnds = betainv([0.025 0.975],bp(1),bp(2));
display(['DRC local transmission relative to a seeding event: ' num2str(mle,'%3.2f') ' (95% CrI:' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')']);

bp = betafit(binopdf(0,2,R(:)));
mle = (bp(1)-1)./(sum(bp)-2);
bnds = betainv([0.025 0.975],bp(1),bp(2));
display(['Probability two cases in DRC only a seeding event: ' num2str(mle,'%3.2f') ' (95% CrI:' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')']);


bp = betafit(binopdf(0,Cases_Uganda_Reported,R_Uganda(:)));
if(bp(1)>1)
    mle = (bp(1)-1)./(sum(bp)-2);
    bnds = betainv([0.025 0.975],bp(1),bp(2));
else
    pu=binopdf(0,Cases_Uganda_Reported,R_Uganda(:));
    [mn,sn] = normfit(log(pu./(1-pu)));
    mle = 1./(1+exp(-mn));
    lb = norminv(0.025,mn,sn);
    ub = norminv(0.975,mn,sn);
    lb=1./(1+exp(-lb));
    ub=1./(1+exp(-ub));
    bnds =[lb ub];
end
display(['Probability ' num2str(Cases_Uganda_Reported) ' cases in Uganda only a seeding event: ' num2str(mle,'%3.2f') ' (95% CrI:' num2str(bnds(1),'%3.2f') char(8211) num2str(bnds(2),'%3.2f') ')']);

[N,edges] = histcounts(Cases_Uganda_Local(:),-0.5:15.5);
case_f = edges(1:end-1)+0.5;
mle = case_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = case_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]);
display(['Number cases associated with localized transmission in Uganda: ' num2str(mle,'%d') ' (95% CrI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);


[N,edges] = histcounts(test_time(:,1),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated: ' datestr(mle) ' (95% CrI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);


[N,edges] = histcounts(test_time(samp_Generation(:,1)==G_ref,1),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated if in generation ' num2str(G_ref,'%d') ': ' datestr(mle) ' (95% CrI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);


display('=================================================================================')
display('Same generation')
display('=================================================================================')

[N,edges] = histcounts(samp_Indx(:,2),0.5:5.5);
mle = find(N==max(N));
test_bnds = cumsum(N)./sum(N);
bnds = [find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")];
display(['Number of index cases: ' num2str(mle,'%d') ' (95% UI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);

[N,edges] = histcounts(samp_Generation(:,2),-0.5:15.5);
G_f = edges(1:end-1)+0.5;
mle = G_f(find(N==max(N)));
mle_Gen = mle;
test_bnds = cumsum(N)./sum(N);
bnds = G_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]);
display(['Generation of epidemic: ' num2str(mle,'%d') ' (95% UI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);


[N,edges] = histcounts(test_time(:,2),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated: ' datestr(mle) ' (95% UI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);

[N,edges] = histcounts(test_time(samp_Generation(:,2)==mle_Gen,2),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated if in generation ' num2str(mle_Gen,'%d') ': ' datestr(mle) ' (95% UI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);

display('=================================================================================')
display('Combined')
display('=================================================================================')

[N,edges] = histcounts(samp_Indx(:,3),0.5:5.5);
mle = find(N==max(N));
test_bnds = cumsum(N)./sum(N);
bnds = [find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")];
display(['Number of index cases: ' num2str(mle,'%d') ' (95% UI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);

[N,edges] = histcounts(samp_Generation(:,3),-0.5:15.5);
G_f = edges(1:end-1)+0.5;
mle = G_f(find(N==max(N)));
mle_Gen = mle;
test_bnds = cumsum(N)./sum(N);
bnds = G_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]);
display(['Generation of epidemic: ' num2str(mle,'%d') ' (95% UI:' num2str(bnds(1),'%d') char(8211) num2str(bnds(2),'%d') ')']);



[N,edges] = histcounts(test_time(:,3),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated: ' datestr(mle) ' (95% UI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);


[N,edges] = histcounts(test_time(samp_Generation(:,3)==mle_Gen,3),-0.5:250);
time_f = edges(1:end-1)+0.5;
mle = datenum('May 13, 2026') - time_f(find(N==max(N)));
test_bnds = cumsum(N)./sum(N);
bnds = datenum('May 13, 2026') -flip(time_f([find(test_bnds>=0.025,1,"first") find(test_bnds>=0.975,1,"first")]));
display(['Date epidemic originated if in generation ' num2str(mle_Gen,'%d') ': ' datestr(mle) ' (95% UI:' datestr(bnds(1)) char(8211) datestr(bnds(2)) ')']);



close all;
figure('units','normalized','outerposition',[0.05 0.05 0.9 1])
subplot('Position',[0.055 0.6 0.94 0.35])


C_Text={'#d0d1e6'
'#a6bddb'
'#67a9cf'
'#3690c0'
'#02818a'
'#016c59'
'#014636'};
C_RGB=flip(hex2rgb(C_Text)); 

num_g =abs(min(samp_Generation(:,1))-max(samp_Generation(:,1)))+1;

CM=interp1([1:size(C_RGB,1)],C_RGB,linspace(1,size(C_RGB,1),num_g));
t_plot=linspace(0,350,10001);

peak_count = zeros(num_g,1);
gen=min(samp_Generation(:,1)):max(samp_Generation(:,1))  ;
for gg=1:num_g
    [N,edges] = histcounts(test_time(samp_Generation(:,1)==gen(gg),1),-0.5:250);
    peak_count(gg)=max(N);
end
[~,S_indx]=sort(peak_count,'descend');
peak_count=peak_count(S_indx);
gen=gen(S_indx);
test_d=datetime('May 13, 2026')-t_plot;
legend_txt=cell(num_g,1);
for gg=1:num_g
    % histogram(test_time(samp_Generation(:,1)==gg,1),-0.5:250); hold on
    gam_plot=gamfit(test_time(samp_Generation(:,1)==gen(gg),1));
    t_mle=(gam_plot(1)-1).*gam_plot(2);
    [N,edges] = histcounts(test_time(samp_Generation(:,1)==gen(gg),1),-0.5:250);
    plot(test_d,max(N).*gampdf(t_plot,gam_plot(1),gam_plot(2))./gampdf(t_mle,gam_plot(1),gam_plot(2)),'color',CM(gg,:),'LineWidth',2); hold on
    legend_txt{gg}=['Generation ' num2str(gen(gg),'%d')];
end

xlim([datetime('October 15, 2025') datetime('May 13, 2026')])
xlabel('Origin date of epidemic','FontSize',16)
ylabel('Likelihood','FontSize',16)
set(gca,'LineWidth',2,'tickdir','out','XTick',flip([datetime('May 13, 2026'):-7:datetime('October 1, 2025')]),'FontSize',12)
box off;
theme('light');
legend(legend_txt,'FontSize',14,'Location','northwest');
legend boxoff;
text(-0.05,1,'A','FontSize',28,'Units','normalized');

subplot('Position',[0.055 0.12 0.94 0.35])
% C_Text={'#bfd3e6'
% '#9ebcda'
% '#8c96c6'
% '#8c6bb1'
% '#88419d'
% '#810f7c'
% '#4d004b'};
C_RGB=flip(hex2rgb(C_Text)); 

num_g =abs(min(samp_Generation(:,2))-max(samp_Generation(:,2)))+1;

CM=interp1([1:size(C_RGB,1)],C_RGB,linspace(1,size(C_RGB,1),num_g));
t_plot=linspace(0,350,10001);

peak_count = zeros(num_g,1);
gen=min(samp_Generation(:,2)):max(samp_Generation(:,2))  ;
for gg=1:num_g
    [N,edges] = histcounts(test_time(samp_Generation(:,2)==gen(gg),2),-0.5:250);
    peak_count(gg)=max(N);
end
[~,S_indx]=sort(peak_count,'descend');
peak_count=peak_count(S_indx);
gen=gen(S_indx);
test_d=datetime('May 13, 2026')-t_plot;
legend_txt=cell(num_g,1);
for gg=1:num_g
    % histogram(test_time(samp_Generation(:,1)==gg,1),-0.5:250); hold on
    gam_plot=gamfit(test_time(samp_Generation(:,2)==gen(gg),2));
    t_mle=(gam_plot(1)-1).*gam_plot(2);
    [N,edges] = histcounts(test_time(samp_Generation(:,2)==gen(gg),2),-0.5:250);
    plot(test_d,max(N).*gampdf(t_plot,gam_plot(1),gam_plot(2))./gampdf(t_mle,gam_plot(1),gam_plot(2)),'color',CM(gg,:),'LineWidth',2); hold on
    legend_txt{gg}=['Generation ' num2str(gen(gg),'%d')];
end

xlim([datetime('October 15, 2025') datetime('May 13, 2026')])
xlabel('Origin date of epidemic','FontSize',16)
ylabel('Likelihood','FontSize',16)
set(gca,'LineWidth',2,'tickdir','out','XTick',flip([datetime('May 13, 2026'):-7:datetime('October 1, 2025')]),'FontSize',12)
box off;
theme('light');
legend(legend_txt,'FontSize',14,'Location','northwest');
legend boxoff;
text(-0.05,1,'B','FontSize',28,'Units','normalized');
print(gcf, 'Supplementary_Figure_7.png','-dpng','-r300');