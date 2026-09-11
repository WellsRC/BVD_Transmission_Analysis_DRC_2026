function gt=Generation_Time(a_gen,b_gen,dist_type,n_samp)
if(strcmp(dist_type,'Gamma'))
    gt=gamrnd(a_gen,b_gen,n_samp,1);
elseif(strcmp(dist_type,'Weibull'))
    gt=wblrnd(a_gen,b_gen,n_samp,1);
end
end