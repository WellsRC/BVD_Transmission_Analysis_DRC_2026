function y = Recursive_Summation_Asymptomatic_Chain(G,current_G,pA,k,p,xf,x_past)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Inputs
% G - Number of generations
% current_G - Current generation
% p_A - asymptomatic proportion
% k - dispersion parameter for the negative binomial
% p - probability parameter for the negative binomial
% x_past - valus of x in past genrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
x_current=[1:nbininv(1-10^(-4),xf*k,p)];
if(current_G==G)
    v = nbinpdf(x_current,xf*k,p).*x_current.*(1-pA).*pA.^(sum(x_past)+x_current);
else
    v=zeros(size(x_current));
    for jj=1:length(v)
        v(jj)=nbinpdf(x_current(jj),xf*k,p).*Recursive_Summation_Asymptomatic_Chain(G,current_G+1,pA,k,p,x_current(jj),[x_past x_current(jj)]);
    end
end

y=sum(v);

end

