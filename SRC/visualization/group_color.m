function c = group_color(genotype, drug_idx, n_drugs, config)
    % GROUP_COLOR - Colour for a genotype x drug group
    %
    % Genotype colour from config.groups; later drugs (drug_idx > 1) are
    % progressively darker shades of it so the genotype stays readable across
    % drug conditions. Near-black bases cannot darken, so they lighten instead.
    idx = find(strcmp({config.groups.name}, genotype), 1);
    if isempty(idx), base = [0.4 0.4 0.4]; else, base = config.groups(idx).color; end
    f = 0.45 * (drug_idx - 1) / max(n_drugs - 1, 1);
    if max(base) < 0.25
        c = base + (1 - base) * f;
    else
        c = base * (1 - f);
    end
end
