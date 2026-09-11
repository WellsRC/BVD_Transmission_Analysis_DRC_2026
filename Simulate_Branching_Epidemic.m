function [Incidence] = Simulate_Branching_Epidemic(k,p,G)

NS = length(k);
Incidence = zeros(NS,G+1);

for ss=1:NS
    Incidence(ss,1)=1;
    for gg=2:G+1
        if(Incidence(ss,gg-1)>0)
            Incidence(ss,gg) = nbinrnd(k(ss).*Incidence(ss,gg-1),p(ss));
        end
    end
end

end
