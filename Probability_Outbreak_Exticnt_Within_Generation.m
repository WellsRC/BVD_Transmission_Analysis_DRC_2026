function [p_extinct,Gen]=Probability_Outbreak_Exticnt_Within_Generation(a_gen,b_gen,r_growth,k,max_gen)

R0=Compute_Reproduction_Number(a_gen,b_gen,r_growth,[]);
p=k./(R0+k);

Gen=[0:max_gen];
p_extinct=zeros(max_gen+1,1);
p_extinct(1)=nbinpdf(0,k,p);
x=[0:10^4];
for jj=2:(max_gen+1)
    p_extinct(jj)=sum(nbinpdf(x,k,p).*p_extinct(jj-1).^x);
end

end