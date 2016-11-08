function Vprinc = studyStochMod(systemName,caseName)
% Study the stoch model of a system
nbPC = 10;
verbose = 0;

stochMod = getSystemStochModel(systemName,caseName,1);
mpc = openCase(systemName);
nzLoads = find(mpc.bus(:,3) ~= 0);
% nbStoch = length(obj.muNorm);
idxStoch = nzLoads;
[V,D] = eig(stochMod.sigNorm);
[eig_sort,idx_sort] = sort(diag(D),'descend');
idxPC = idx_sort(1:min(nbPC,length(idx_sort)));
Vprinc = V(:,idxPC);
idxNZinV = find(sum(Vprinc,2) ~= 0);

if verbose
    % Legend
    legendTxt = cell(nbPC,1);
    for i = 1:nbPC
        legendTxt{i} = sprintf('For eigenvalue %d',i);
    end
    
    figure
    plot(eig_sort,'*');
    title('Eigenvalues in descending order');
    
    figure
    plot(nzLoads,Vprinc.','LineWidth',2);
    legend(legendTxt);
    title('Components in eigenvectors');
    
    figure
    hold on
    bar(nzLoads,stochMod.muNorm);
    bar(nzLoads(idxNZinV),stochMod.muNorm(idxNZinV),'g');
    title('Mean value of load');
    
    figure
    surf(stochMod.sigNorm(idxNZinV,idxNZinV));
    title('Covariance matrix of the five PC');
end