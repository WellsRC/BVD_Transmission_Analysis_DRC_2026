clear;
clc;
Max_Num_Gen =[15 15 15 15 15];
rng(20260720)

[k_mle_base,p_mle,cfr_mle,k,p,cfr,R0_mle,R0] = BVD_Parameters(1);
pA_mle = 1-cfr_mle;
Lt_Lagged = zeros(5,max(Max_Num_Gen)+1); % likelihood using a lag for the deaths; ie deaths based on previous generations only (i.e. if at generation G looking at 0:G-1)
Lt_Generation = zeros(5,max(Max_Num_Gen)+1); % likelihood is based on deaths and incidence in the same generation


N_Samp=10^4;

for Num_Index_Cases=1:5
    Num_Gen= Max_Num_Gen(Num_Index_Cases);
    % Run the branching process for the simualted epidemic
    for ss=1:Num_Index_Cases
        if(ss==1)
            [Incidence] = Simulate_Branching_Epidemic(k_mle_base.*ones(N_Samp,1),p_mle.*ones(N_Samp,1),Num_Gen);
        else            
            Incidence = Incidence+Simulate_Branching_Epidemic(k_mle_base.*ones(N_Samp,1),p_mle.*ones(N_Samp,1),Num_Gen);
        end
    end
    % Adjust the dispersion coefficient based on the number of index cases
    k_mle = k_mle_base.*Num_Index_Cases;

    % Calculate the number of deaths
    Death = zeros(size(Incidence)) ;    
    for ss=1:N_Samp
        for gg=1:Num_Gen+1
            if(Incidence(ss, gg).*cfr_mle.*(1-cfr_mle)>=10^6)
                Death(ss,gg) =round(normrnd(Incidence(ss, gg).*cfr_mle,sqrt(Incidence(ss, gg).*cfr_mle.*(1-cfr_mle))));
            else
                Death(ss,gg) = sum(rand(Incidence(ss, gg),1)<= cfr_mle);
            end
        end
    end

    % Calcualte the likelihood/pdf of detecting S_mle deaths
    S_mle=80;
    mle_Detect_Chain=Sample_Approximation_Detection_Chain(Num_Gen,S_mle,pA_mle,k_mle,p_mle,Num_Index_Cases);    
    pdf_Detect = [mle_Detect_Chain(1); diff(mle_Detect_Chain)];


    Gen =[0:Num_Gen];    
    Case_Data=246;
    Death_Data =S_mle;

    Death_Cumulative = cumsum(Death,2);
    Incidence_Cumulative = cumsum(Incidence,2);

    L_Lag = zeros(max(Max_Num_Gen)+1,3); % Likelihood for deaths based on the prior generation
    L_G = zeros(max(Max_Num_Gen)+1,3); % Likelihood for deaths and incdience based on the same generation
    for gg=2:Num_Gen+1        
        nb_case = nbinfit(Incidence_Cumulative(Incidence(:,gg-1)>0,gg));
        L_Lag(gg,1) = nbinpdf(Case_Data,nb_case(1),nb_case(2));
        L_G(gg,1) = nbinpdf(Case_Data,nb_case(1),nb_case(2));
        if(gg>2)
            nb_inc = nbinfit(Incidence(:,gg-1));
            L_Lag(gg,1) = nbinpdf(Case_Data,nb_case(1),nb_case(2))./(1-nbinpdf(0,nb_inc(1),nb_inc(2)));
            L_G(gg,1) = nbinpdf(Case_Data,nb_case(1),nb_case(2))./(1-nbinpdf(0,nb_inc(1),nb_inc(2)));
        end
        if(gg>3)
            nb_inc = nbinfit(Incidence(:,gg-2));
            nb_death = nbinfit(Death_Cumulative(Incidence(:,gg-2)>0,gg-1)); % Need to look at the numebr of deaths coming from the previous generation as deaths in current generation have not occurred yet
            L_Lag(gg,2) = nbinpdf(Death_Data,nb_death(1),nb_death(2))./(1-nbinpdf(0,nb_inc(1),nb_inc(2)));
        else
            nb_death = nbinfit(Death_Cumulative(:,gg-1)); % Need to look at the numebr of deaths coming from the previous generation as deaths in current generation have not occurred yet
            if(~isinf(nb_death(1)))
                L_Lag(gg,2) = nbinpdf(Death_Data,nb_death(1),nb_death(2));
            else
                nb_death =poissfit(Death_Cumulative(:,gg-1));
                L_Lag(gg,2) = poisspdf(Death_Data,nb_death);
            end
        end
        if(gg>2)
            nb_inc = nbinfit(Incidence(:,gg-1));
            nb_death = nbinfit(Death_Cumulative(Incidence(:,gg-1)>0,gg))./(1-nbinpdf(0,nb_inc(1),nb_inc(2))); % Need to look at the numebr of deaths coming from the previous generation as deaths in current generation have not occurred yet
        else
            nb_death = nbinfit(Death_Cumulative(Incidence(:,gg-1)>0,gg));
        end
        L_G(gg,2) = nbinpdf(Death_Data,nb_death(1),nb_death(2));

        L_Lag(gg,3)=pdf_Detect(gg-1); %Detection is based on deaths, so need to look at the prior generation
        L_G(gg,3)=pdf_Detect(gg);

    end
    Lt_Lagged(Num_Index_Cases,:)=prod(L_Lag,2);
    Lt_Generation(Num_Index_Cases,:)=prod(L_G,2);
end

save('Epidemic_Analysis_Likelihood');