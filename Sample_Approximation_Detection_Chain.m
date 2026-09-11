function P_detect = Sample_Approximation_Detection_Chain(G,detect_S,pA,k,p,I0)

P_detect = zeros(G+1,1);
P_detect(1)=1-binocdf(detect_S-1,1,1-pA);
if(G>0)
    N_Samp=10^4;
    Samp_Gen_Inf = zeros(G,N_Samp);
    for gg=1:G
        if(gg==1)
            x_past=ones(1,N_Samp);
        else
            x_past = Samp_Gen_Inf(gg-1,:);
        end
        test_samp=nbinrnd(x_past.*k,p,1,N_Samp);
        while(sum(test_samp==0)>0)
            indx_rs = test_samp==0;
            test_samp(indx_rs)=nbinrnd(x_past(indx_rs).*k,p,1,sum(indx_rs));
        end
        Samp_Gen_Inf(gg,:)=test_samp;
    end
    Num_I = I0+cumsum(Samp_Gen_Inf(1:gg,:),1);
    dX=Num_I-(detect_S-1);
    y=ones(size(dX));
    y(dX>0) = betainc(1-(1-pA), dX(dX>0), (detect_S-1)+1);
    P_detect(2:end) = mean(1-y,2);
end

end
