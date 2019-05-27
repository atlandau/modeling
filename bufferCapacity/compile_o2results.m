hPath = '/Users/landauland/Documents/Research/SabatiniLab/modeling/bufferCapacity';
rPath = fullfile(hPath,'results_run190523');
resultDir = dir(rPath);
resultName = {resultDir(:).name};
resultName = resultName(contains(resultName,'results_'));

% I just know these...
NK = 9; % kon
ND = 7; % kd
NA = 7; % amplitude

% Timing Parameters - forgot to save these in the output parameter struct
dt = 1e-8;
ds = 1000;
T = 0.3;
tvec = 0:dt*ds:T;

p = cell(NK,ND,NA); % parameter structure from each run
y = cell(NK,ND,NA); % data from each run
dk = cell(NK,ND,NA); % dynamic kappa, computed here

for nk = 1:NK
    fprintf('Kon: %d/%d...\n',nk,NK);
    for nd = 1:ND
        for na = 1:NA
            cName = sprintf('results_Kon%d_Kd%d_Amp%d.mat',nk,nd,na);
            cResult = load(fullfile(rPath,cName));
            p{nk,nd,na} = cResult.p;
            y{nk,nd,na} = cResult.y;
            dk{nk,nd,na} = cResult.y(:,2) ./ cResult.y(:,1);
            clear cResult
        end
    end
end


            


