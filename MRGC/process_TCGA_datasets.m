function [U, idx_eg, idx_rc] = process_TCGA_datasets(X, par)
    U = MRGC(X, par.num_bases, par.alpha, par.beta, par.num_iters); % X is a cell of matrix, each matrix is n*m, n is the number of samples, m is the number of features, each row is a data point.
    NUMC = 2 : 15;
    [K1, ~, K12, ~] = Estimate_Number_of_Clusters_given_graph(U, NUMC);
    fprintf('The number of clusters estimated by eigen gap : %d\n', K1);
    idx_eg = computeLabels(U, K1);
    fprintf('The number of clusters estimated by rotation cost : %d\n', K12);
    idx_rc = computeLabels(U, K12);
end


function idx = computeLabels(T, k)
    Z = (abs(T)+ abs(T')) / 2;
    idx = clu_ncut(Z, k);
end
