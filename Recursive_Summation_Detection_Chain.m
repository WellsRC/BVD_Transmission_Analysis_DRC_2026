function y = Recursive_Summation_Detection_Chain(G,I_max,current_G,S,pA,k,p,xG_past,x_sum_past)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Inputs
% G - Number of generations
% current_G - Current generation
% S - d
% p_A - asymptomatic proportion
% k - dispersion parameter for the negative binomial
% p - probability parameter for the negative binomial
% x_past - valus of x in past genrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if(G==0)
    v=(1-binocdf(S-1,1,1-pA));
elseif(current_G==1 && G==1)
    x_current=[1:I_max];
    v=nbinpdf(0,k,p).*(1-binocdf(S-1,1,1-pA))+(1-nbincdf(I_max,k,p))+nbinpdf(x_current,k,p).*(1-binocdf(S-1,1+x_current,1-pA));
elseif(current_G==G-1)
    if(I_max-sum(x_sum_past)>=1)
        x_current=[1:I_max-sum(x_sum_past)];
        temp_next_Gen=zeros(size(x_current));
        for jj=1:length(x_current)
            x_temp=[0:I_max-sum(x_sum_past)-x_current(jj)];
            if(~isempty(x_temp))
                temp_next_Gen(jj)=sum(nbinpdf(x_temp,x_current(jj).*k,p).*(1-binocdf(S-1,1+sum(x_sum_past)+x_current(jj)+x_temp,1-pA)));
            end
        end
        v = nbinpdf(x_current,xG_past.*k,p).*(temp_next_Gen-nbincdf(I_max-sum(x_sum_past)-x_current,x_current.*k,p));
        if(G-1==1)
            temp_C = nbinpdf(0,k,p).*(1-binocdf(S-1,1,1-pA))+(1-nbincdf(0,k,p));
            v=sum(v)+temp_C;
        end
    else
        v=0;
    end
elseif(current_G==1)
    temp_C = nbinpdf(0,k,p).*(1-binocdf(S-1,1,1-pA))+(1-nbincdf(0,k,p));
    x_current=[1:I_max];
    v=zeros(size(x_current));
    for jj=1:length(x_current)
        v(jj) = nbinpdf(x_current(jj),k,p).*Recursive_Summation_Detection_Chain(G,I_max,current_G+1,S,pA,k,p,x_current(jj),[x_sum_past x_current(jj)]);
    end
    v = sum(v)+temp_C;
else
    x_current=[1:I_max-sum(x_sum_past)];
    if(I_max-sum(x_sum_past)>=1)
        v=zeros(size(x_current));
        parfor jj=1:length(x_current)
            v(jj) = nbinpdf(x_current(jj),xG_past.*k,p).*Recursive_Summation_Detection_Chain(G,I_max,current_G+1,S,pA,k,p,x_current(jj),[x_sum_past x_current(jj)]);
        end
    else
        v=0;
    end
end


y=sum(v);

end

