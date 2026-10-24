function compat = check_pyenv(py2ok)

if nargin < 1
    py2ok = false;
end

rel.R2017b = {'2.7', '3.4', '3.5', '3.6'};
rel.R2018a = {'2.7', '3.5', '3.6'};
rel.R2018b = {'2.7', '3.5', '3.6'};
rel.R2019a = {'2.7', '3.5', '3.6', '3.7'};
rel.R2019b = {'2.7', '3.6', '3.7'};
rel.R2020a = {'2.7', '3.6', '3.7'};
rel.R2020b = {'2.7', '3.6', '3.7', '3.8'};
rel.R2021a = {'2.7', '3.7', '3.8'};
rel.R2021b = {'2.7', '3.7', '3.8', '3.9'};
rel.R2022a = {'2.7', '3.8', '3.9'};
rel.R2022b = {'2.7', '3.8', '3.9', '3.10'};
rel.R2023a = {'3.8', '3.9', '3.10'};

v = ver('MATLAB');
mlver = v.Release(2:7);
mlyear = str2double(mlver(2:5));
mlab = mlver(6);

if mlyear < 2019
    pyver = pyversion.Rersion;
elseif mlyear == 2019
    if strcmp(mlab, 'a')
        pyver = pyversion.Rersion;
    else
        pyver = pyenv().Version;
    end
else
    pyver = pyenv().Version;
end

if any(ismember(rel.(mlver), pyver))
    compat = true;
else
    compat = false;
end

if pyver == '2.7' && py2ok == false
    compat = false;
end

