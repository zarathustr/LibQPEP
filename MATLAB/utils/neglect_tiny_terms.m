% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Unified Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description 
%           in A Globally Optimal Framework.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com


function new_eq = neglect_tiny_terms(eqs, level)
thresh = 10.^(-level);
eqs = vpa(eqs, 64);
for i = 1 : length(eqs)
    terms = children(eqs(i));
    sum = 0;
    for j = 1 : length(terms)
        term = terms(j);
        facts = children(term);
        fact = facts(end);
        str = char(fact);
        str_val = double(str(1));
        if(~((str_val >= double('A') && str_val <= double('Z')) || (str_val >= double('a') && str_val <= double('z'))))
            if(abs(fact) <= thresh)
                term = 0;
            end
        end
        sum = sum + term;
    end
    new_eq(i) = sum;
end
end