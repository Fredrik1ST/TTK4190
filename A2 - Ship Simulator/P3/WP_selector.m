function [xk1, yk1, xk, yk, last, k] = WP_selector(xN, yE, WP)
% WP_SELECTOR: selects active waypoint segment using input WP
% Inputs:
%   xN, yE : current North/East position [m]
%   WP     : waypoints as 2×N (rows: North; East) or N×2
% Outputs:
%   (xk,yk)->(xk1,yk1) : current active segment
%   last   : true if final waypoint reached (within Rsw)
%   k      : current segment index (between 1 and N-1)
%
% Notes:
% - No file I/O. Uses WP passed as argument.
% - Persists the segment index k across calls. If WP size changes, k resets.

    % --- Validate/reshape WP to 2×N ---
    if size(WP,1) ~= 2
        if size(WP,2) == 2
            WP = WP.'; % make it 2×N if N×2
        else
            error('Waypoints must be 2×N (rows: North; East) or N×2.');
        end
    end

    % --- Persistent state: current segment index and switching radius ---
    persistent k_persist nWP_prev Rsw
    if isempty(k_persist) || isempty(nWP_prev) || nWP_prev ~= size(WP,2)
        k_persist = 1;              % reset if first call or WP changed
        nWP_prev  = size(WP,2);
        Rsw       = 500;            % switching radius [m] — adjust if needed
    end

    % --- Bounds / basic checks ---
    nWP = size(WP,2);
    if nWP < 2
        error('Need at least two waypoints.');
    end
    if k_persist >= nWP
        k_persist = nWP - 1;
    end

    % --- Current segment endpoints ---
    xk  = WP(1,k_persist);   yk  = WP(2,k_persist);
    xk1 = WP(1,k_persist+1); yk1 = WP(2,k_persist+1);

    % --- Geometry along the segment ---
    dx = xk1 - xk;  dy = yk1 - yk;
    L  = hypot(dx, dy) + 1e-9; % segment length (avoid division by zero)

    % Along-track projection of current position onto segment direction
    s = ((xN - xk)*dx + (yE - yk)*dy) / L;

    % --- Switch to next segment when close to the end of the current one ---
    if s >= max(L - Rsw, 0) && k_persist < nWP-1
        k_persist = k_persist + 1;
        xk  = WP(1,k_persist);   yk  = WP(2,k_persist);
        xk1 = WP(1,k_persist+1); yk1 = WP(2,k_persist+1);
        dx = xk1 - xk; dy = yk1 - yk; L = hypot(dx,dy) + 1e-9;
    end

    % --- Final waypoint reached? ---
    last = (k_persist == nWP-1) && (hypot(xN - xk1, yE - yk1) < Rsw);

    % Expose current segment index
    k = k_persist;
end
