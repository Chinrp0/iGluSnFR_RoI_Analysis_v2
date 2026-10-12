function [geno_list, drug_list] = genotype_drug_order(genotypes, drugs, config)
    % GENOTYPE_DRUG_ORDER - Genotypes and drugs present, in plotting order
    %
    % Genotypes follow config.groups order; drugs follow config.drugs.order,
    % then any drug not listed there (alphabetical). Shared by the 1AP and
    % spontaneous batch drivers so both form genotype x drug groups the same way.

    cfg_geno  = {config.groups.name};
    geno_list = cfg_geno(ismember(cfg_geno, genotypes));
    cfg_drug  = config.drugs.order;
    drug_list = cfg_drug(ismember(cfg_drug, drugs));
    extra     = setdiff(unique(drugs), cfg_drug);
    drug_list = [drug_list(:)', extra(:)'];
end
