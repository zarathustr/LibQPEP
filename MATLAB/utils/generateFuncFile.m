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


function generateFuncFile(obj, filename, vars)
if(verLessThan('matlab', '8.5.0'))
    matlabFunction(obj, 'File', filename, 'Vars', vars);
else
    matlabFunction(obj, 'File', filename, 'Optimize', true, 'Vars', vars);
end

fp = fopen(filename, 'r');
counter = 0;
while(true)
    ll = fgetl(fp);
    if(~isempty(ll))
        if(isempty(strfind(ll, '%')))
            counter = counter + 1;
            lls{counter} = ll;
        end
    end
    
    if(ll == -1)
        break;
    end
end
fclose(fp);

fp = fopen(filename, 'w+');
for i = 1 : length(lls) - 1
    fprintf(fp, '%s\n', lls{i});
end
fclose(fp);
end