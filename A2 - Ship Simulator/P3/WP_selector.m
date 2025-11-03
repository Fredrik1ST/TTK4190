function [xk1, yk1, xk, yk, last, k] = WP_selector(xN, yE, WP, Rsw, Rstop)
% WP_SELECTOR: Selects the active waypoint segment and detects when the
% final waypoint has been reached.
%
% Inputs:
%   xN, yE : current North/East position [m]
%   WP     : waypoints as 2×N (rows: North; East) or N×2
%   Rsw    : switching radius [m] (optional, default 300)
%   Rstop  : stop radius for final waypoint [m] (optional, default 50)
%
% Outputs:
%   (xk,yk)->(xk1,yk1) : current active segment
%   last   : true if final waypoint reached (within Rstop)
%   k      : current segment index (1..N-1)
%
% Notes:
% - Fully general for any number of waypoints N ≥ 2
% - Uses persistent k between calls
% - Automatically resets if WP array changes

    % --- Defaults ---
    if nargin < 4 || isempty(Rsw),   Rsw = 300; end
    if nargin < 5 || isempty(Rstop), Rstop = 50;  end

    % --- Validate/reshape WP to 2×N ---
    if size(WP,1) ~= 2
        if size(WP,2) == 2
            WP = WP.';  % transpose if N×2
        else
            error('Waypoints must be 2×N or N×2.');
        end
    end
    nWP = size(WP,2);
    if nWP < 2
        error('Need at least two waypoints.');
    end

    % --- Persistent: current segment index ---
    persistent k_persist nWP_prev
    if isempty(k_persist) || isempty(nWP_prev) || nWP_prev ~= nWP
        % reset if first call or WP changed
        k_persist = 1;
        nWP_prev  = nWP;
    end
    if k_persist >= nWP
        k_persist = nWP - 1;
    end

    % --- Helper to get current segment info ---
    function [xk, yk, xk1, yk1, dx, dy, L] = segment(kidx)
        xk  = WP(1,kidx);    yk  = WP(2,kidx);
        xk1 = WP(1,kidx+1);  yk1 = WP(2,kidx+1);
        dx = xk1 - xk; dy = yk1 - yk;
        L  = hypot(dx, dy) + 1e-12;
    end

    % --- Progress to next segment if near end ---
    advanced = true;
    while advanced
        advanced = false;
        [xk, yk, xk1, yk1, dx, dy, L] = segment(k_persist);

        % Distance metrics
        s = ((xN - xk)*dx + (yE - yk)*dy) / L;  % projection (m)
        dist_to_end = hypot(xN - xk1, yE - yk1);

        % Move to next segment if within switch radius
        if k_persist < nWP-1 && (dist_to_end < Rsw || s >= L - Rsw)
            k_persist = k_persist + 1;
            advanced = true;
        end
    end

    % --- Output current segment ---
    [xk, yk, xk1, yk1, ~, ~, ~] = segment(k_persist);

    % --- Check for final waypoint ---
    dist_final = hypot(xN - WP(1,end), yE - WP(2,end));
    last = (k_persist == nWP-1) && (dist_final <= Rstop);

    % --- Return current segment index ---
    k = k_persist;
end
