function eta_G=Probability_Asymptomatic_Index_Symptomatic_Generation(p,k,p_A,eta_Total)
    eta_G = p_A-eta_Total-p_A.*(p.^k).*(1-(p_A-eta_Total).*(1-p)).^(-k);

    eta_G=max(eta_G,0);
end