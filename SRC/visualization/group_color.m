function c = group_color(genotype, drug_idx, n_drugs, config)
    % GROUP_COLOR - Colour for a genotype x drug group
    %
    % Genotype colour from config.groups; later drugs (drug_idx > 1) are
    % progressively darker shades of it so the genotype stays readable across
    % drug conditions.
    idx = find(strcmp({config.groups.name}, genotype), 1);
    if isempty(idx), base = [0.4 0.4 0.4]; else, base = config.groups(idx).color; end
    shade = 1 - 0.45 * (drug_idx - 1) / max(n_drugs - 1, 1);
    c = base * shade;
end
