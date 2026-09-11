function [a_gen,b_gen,dist_type]=Parameter_Generation_Time(Disease)

if(strcmp(Disease,'Influenza'))
%https://www.nejm.org/doi/suppl/10.1056/NEJMoa0906089/suppl_file/nejm_lessler_2628sa1.pdf
    a_gen=2.36;
    b_gen=3.18;
    dist_type='Weibull';
elseif(strcmp(Disease,'Monkeypox'))
    % https://pmc.ncbi.nlm.nih.gov/articles/instance/9514338/bin/22-1126-Techapp-s1.pdf
    a_gen=4.85;
    b_gen=2.57;
    dist_type='Gamma';
elseif(strcmp(Disease,'Ebola'))
%https://pmc.ncbi.nlm.nih.gov/articles/PMC4207625/
    a_gen=13.60;
    b_gen=2.59;
    dist_type='Weibull';
elseif(strcmp(Disease,'COVID'))
    % https://www.sciencedirect.com/science/article/pii/S2468042720300634?via%3Dihub
    
    a_gen=1.813;
    b_gen=2.199;
    dist_type='Gamma';
end