function R0=Compute_Reproduction_Number(a_gen,b_gen,dist_type,r_growth,tau)
if(strcmp(dist_type,'Gamma'))
    if(isempty(tau))
        M=integral(@(t)exp(-r_growth.*t).*gampdf(t,a_gen,b_gen),0,inf,'AbsTol',10^(-14),'RelTol',10^(-14));
        if(isnan(M) || isinf(M))
            MF=4;
            while(isnan(M) || isinf(M))
                M=integral(@(t)exp(-r_growth.*t).*gampdf(t,a_gen,b_gen),0,10^MF);
                MF=MF-0.1;
            end
        end
        R0=1./M;
    else
        M=integral(@(t)exp(-r_growth.*t).*gampdf(t,a_gen,b_gen),0,inf,'AbsTol',10^(-14),'RelTol',10^(-14));
        if(isnan(M) || isinf(M))
            MF=4;
            while(isnan(M) || isinf(M))
                M=integral(@(t)exp(-r_growth.*t).*gampdf(t,a_gen,b_gen),0,10^MF);
                MF=MF-0.1;
            end
        end
        M_tau=gamcdf(tau,a_gen,b_gen);
        R0=(1./M).*M_tau;
    end
elseif(strcmp(dist_type,'Weibull'))
    if(isempty(tau))
        M=integral(@(t)exp(-r_growth.*t).*wblpdf(t,a_gen,b_gen),0,inf,'AbsTol',10^(-14),'RelTol',10^(-14));
        if(isnan(M) || isinf(M))
            MF=4;
            while(isnan(M) || isinf(M))
                M=integral(@(t)exp(-r_growth.*t).*wblpdf(t,a_gen,b_gen),0,10^MF,'AbsTol',10^(-14),'RelTol',10^(-14));
                MF=MF-0.1;
            end
        end
        R0=1./M;
    else
        M=integral(@(t)exp(-r_growth.*t).*wblpdf(t,a_gen,b_gen),0,inf,'AbsTol',10^(-14),'RelTol',10^(-14));
        if(isnan(M) || isinf(M))
            MF=4;
            while(isnan(M) || isinf(M))
                M=integral(@(t)exp(-r_growth.*t).*wblpdf(t,a_gen,b_gen),0,10^MF,'AbsTol',10^(-14),'RelTol',10^(-14));
                MF=MF-0.1;
            end
        end
        M_tau=wblcdf(tau,a_gen,b_gen);
        R0=(1./M).*M_tau;
    end
end
end