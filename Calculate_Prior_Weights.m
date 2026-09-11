function weight_NEW = Calculate_Prior_Weights(weight_OLD,q_Parameter_MAT,L_Prior)
if(sum(isnan(weight_OLD))>0)
       weight_NEW =1;       
else
       weight_NEW = exp(L_Prior)./(weight_OLD*prod(q_Parameter_MAT,2));
end
end