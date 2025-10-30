function [xk1, yk1, xk, yk, last] = WP_selector(xN, yE, WP)
% WP_SELECTOR: selects active waypoint segment using WP.mat
% Inputs:  xN,yE = current North/East [m]
% Outputs: (xk,yk)->(xk1,yk1) current segment, last=true if final reached

    persistent W k Rsw
    if isempty(W)
        S = load('WP.mat');
        if isfield(S,'W')
            W = S.W;
        elseif isfield(S,'WP')
            W = S.WP;
        else
            error('WP.mat must contain variable named W or WP (2×N).');
        end
        if size(W,1) ~= 2
            if size(W,2) == 2
                W = W.';    % make it 2×N if N×2
            else
                error('Waypoints must be 2×N (rows: North; East).');
            end
        end
        k   = 1;
        Rsw = 500;    % switching radius [m] — juster ved behov
    end

    nWP = size(W,2);
    if nWP < 2, error('Need at least two waypoints.'); end
    if k >= nWP, k = nWP-1; end

    xk  = W(1,k);   yk  = W(2,k);
    xk1 = W(1,k+1); yk1 = W(2,k+1);

    dx = xk1-xk;  dy = yk1-yk;
    L  = hypot(dx,dy) + 1e-9;

    % along-track projection
    s = ((xN - xk)*dx + (yE - yk)*dy) / L;

    % switch near the end of current segment
    if s >= max(L - Rsw, 0) && k < nWP-1
        k   = k + 1;
        xk  = W(1,k);   yk  = W(2,k);
        xk1 = W(1,k+1); yk1 = W(2,k+1);
    end

    last = (k == nWP-1) && (hypot(xN - xk1, yE - yk1) < Rsw);
end
