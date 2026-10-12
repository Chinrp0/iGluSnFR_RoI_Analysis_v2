function [tau_fast, tau_slow, A_fast, A_slow, C, r2] = fit_decay_double(t, y)
    % FIT_DECAY_DOUBLE - Double-exponential decay fit
    %   y(t) ~ A_fast*exp(-t/tau_fast) + A_slow*exp(-t/tau_slow) + C
    % Curve Fitting Toolbox nonlinear least squares, both amplitudes >= 0 and
    % both time constants > 0. Returned with the FAST component first
    % (tau_fast < tau_slow) so the slow tail is unambiguous. Use on a high-SNR
    % averaged trace; single ROIs lack the SNR to resolve two components.
    %
    % Inputs:
    %   t - [n x 1] time from the decay onset (ms), starting at 0
    %   y - [n x 1] baseline-subtracted decay (dF/F), decaying toward ~0
    %
    % Outputs:
    %   tau_fast, tau_slow - time constants (ms), fast first
    %   A_fast, A_slow     - component amplitudes (dF/F)
    %   C                  - offset
    %   r2                 - fit R^2

    tau_fast = NaN; tau_slow = NaN; A_fast = NaN; A_slow = NaN; C = NaN; r2 = NaN;

    good = isfinite(t) & isfinite(y);
    t = double(t(good)); y = double(y(good));
    t = t(:); y = y(:);
    if numel(y) < 6
        return;
    end

    span = max(t(end) - t(1), eps);
    A0   = max(max(y) - min(y), eps);
    C0   = min(y);
    % Start the two components an order of magnitude apart (fast ~ a few frames,
    % slow ~ a sizeable fraction of the trace) so the optimiser separates them.
    ft = fittype('Af*exp(-x/tf) + As*exp(-x/ts) + C', 'independent', 'x', ...
        'coefficients', {'Af', 'tf', 'As', 'ts', 'C'});
    opts = fitoptions(ft);
    opts.StartPoint = [0.7*A0, max(span/20, 5), 0.3*A0, max(span/3, 20), C0];
    opts.Lower      = [0, 1e-3, 0, 1e-3, -Inf];   % amplitudes >= 0, taus > 0
    opts.Upper      = [Inf, Inf, Inf, Inf, Inf];
    opts.MaxIter    = 1000;

    try
        [fobj, gof] = fit(t, y, ft, opts);
    catch
        return;                              % failed to converge
    end

    if ~(isfinite(fobj.tf) && isfinite(fobj.ts) && fobj.tf > 0 && fobj.ts > 0)
        return;
    end
    % Order fast first
    if fobj.tf <= fobj.ts
        tau_fast = fobj.tf; A_fast = fobj.Af;
        tau_slow = fobj.ts; A_slow = fobj.As;
    else
        tau_fast = fobj.ts; A_fast = fobj.As;
        tau_slow = fobj.tf; A_slow = fobj.Af;
    end
    C  = fobj.C;
    r2 = gof.rsquare;
end
