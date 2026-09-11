function  [G,C,Asymp,Gen_Time,I_t,Gen_Count]=NB_Branching_Process(a_gen,b_gen,dist_type,k,pA,R0,R0_tau)

Outbreak_Size=10^4;

p=k./(R0+k);
p_tau=k./(R0_tau+k);

G=zeros(Outbreak_Size);
Gen_Time=zeros(Outbreak_Size);
C=NaN.*zeros(Outbreak_Size,3);
Asymp=zeros(Outbreak_Size,1);
I_t=zeros(Outbreak_Size,1);
if(rand(1)<pA)
    I=nbinrnd(k,p);
    C(1,:)=[0 0 0];
    Asymp(1)=1;
else
    I=nbinrnd(k,p_tau);
    C(1,:)=[1 0 0];
end
Gen_Count=0;
if(I>0)
    G(1,1+[1:I])=1;
    f_indx=find(G(1,:)==1);
    gt=Generation_Time(a_gen,b_gen,dist_type,I);
    Gen_Time(1,f_indx)=gt;
    I_t(f_indx)=gt;
    while(~isempty(f_indx) && find(sum(G,1)>0, 1, 'last' )<Outbreak_Size)
        Gen_Count=Gen_Count+1;
        g_indx=[];
        for jj=1:length(f_indx)
            mx=find(sum(G,1)>0, 1, 'last' );
            if(rand(1)<pA)
                I=nbinrnd(k,p);
                C(f_indx(jj),:)=[0 0 0];
                Asymp(1)=1;
            else
                I=nbinrnd(k,p_tau);
                C(f_indx(jj),:)=[1 0 0];
            end
            if(I>0)
                if(mx+I<=Outbreak_Size)
                    G(f_indx(jj),mx+[1:I])=1;
                    g_indx=[g_indx mx+[1:I]];
                    gt=Generation_Time(a_gen,b_gen,dist_type,I);
                    Gen_Time(f_indx(jj),mx+[1:I])=gt;
                    I_t(mx+[1:I])=I_t(f_indx(jj))+gt;
                else
                    G(f_indx(jj),(mx+1):end)=1;
                    gt=Generation_Time(a_gen,b_gen,dist_type,length(G(f_indx(jj),(mx+1):end)));
                    Gen_Time(f_indx(jj),(mx+1):end)=gt;
                    I_t((mx+1):end)=I_t(f_indx(jj))+gt;
                    break;
                end
            end
        end
        f_indx=g_indx;
    end
    fnz=find(sum(G,1)>0, 1, 'last' );
    G=G(1:fnz,1:fnz);
    C=C(1:fnz,:);
    Gen_Time=Gen_Time(1:fnz,1:fnz);
    I_t=I_t(1:fnz);
    C=C(1:fnz,:);
    f_nan=isnan(C(:,1));
    asymp_case=rand(sum(f_nan),1)<pA;
    Ct=zeros(sum(f_nan),3);
    Ct(~asymp_case,:)=repmat([1 0 0],sum(~asymp_case),1);
    C(f_nan,:)=Ct;
    pAt=ones(sum(f_nan),1);
    pAt(~asymp_case)=0; % put as zero since symptomatic
    Asymp(f_nan)=pAt;
    Asymp=Asymp(1:fnz);
else
    G=0;
    I_t=0;
    Gen_Time=0;
    C=C(1,:);
    Asymp=Asymp(1);
end

end

