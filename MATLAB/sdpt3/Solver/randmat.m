%%******************************************************
%% randmat: generate an mxn matrix using matlab's
%%          rand or randn functions using state = k.
%%
%%******************************************************
function v = randmat(m,n,k,randtype)

try
    s = rng;
    rng(k);
    if strcmp(randtype,'n')
        v = randn(m,n);
    elseif strcmp(randtype,'u')
        v = rand(m,n);
    end
    rng(s);
catch
    if strcmp(randtype,'n')
        s = randn('state'); %#ok
        randn('state',k); %#ok
        v = randn(m,n);
        randn('state',s); %#ok
    elseif strcmp(randtype,'u')
        s = rand('state'); %#ok
        rand('state',k); %#ok
        v = rand(m,n);
        rand('state',s); %#ok
    end
end
%%******************************************************