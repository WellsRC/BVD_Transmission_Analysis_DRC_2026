clear;
clc;
% parpool(6);
N_Samp=10^4;

[k_mle,p_mle,p_fatal_mle,p_non_fatal_mle,cfr_mle,k,p,p_fatal,p_non_fatal,cfr,R0_mle,R0,R0_Fatal_mle,R0_Fatal,R0_Non_fatal_mle,R0_Non_fatal] = BVD_Parameters(N_Samp);

Num_Gen=6;
% uv = [0.05 0.1:0.1:0.9 0.95];
theta_G=zeros(Num_Gen+1,1); 
sigma_G_1=zeros(Num_Gen+1,101); 

theta_G(1)=nbinpdf(0,k_mle,p_mle); % Probability extinct w/ index

for gg=2:Num_Gen+1
    theta_G(gg)=p_mle.^k_mle./(1-sum(theta_G(1:(gg-1))).*(1-p_mle)).^k_mle-sum(theta_G(1:(gg-1)));
end

R0v=flip(linspace(1.01,10,101));
uv=1-1./R0v;
p_v=k_mle./(k_mle+R0v);
for uu=1:length(R0v)
    sigma_G_1(:,uu)=Sample_Approximation_Detection_Chain(Num_Gen,1,uv(uu),k_mle,p_v(uu));
end