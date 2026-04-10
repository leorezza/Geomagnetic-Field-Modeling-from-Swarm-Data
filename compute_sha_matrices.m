function [A_r, A_th, A_ph] = compute_sha_matrices(radius, colat, lon, max_deg, varargin)
% COMPUTE_SHA_MATRICES Generates design matrices for spherical harmonic expansion.
% Supports internal ('int') and external ('ext') field sources.

    % Parse input source type
    source_type = 'int';
    if nargin > 4
        source_type = varargin{1};
    end

    % Ensure inputs are consistent column vectors
    sz = max([size(radius); size(colat); size(lon)]);
    num_pts = prod(sz);
    
    if isscalar(radius); radius = radius * ones(sz); end
    if isscalar(colat); colat = colat * ones(sz); end
    if isscalar(lon); lon = lon * ones(sz); end
    
    r = radius(:);
    th = colat(:);
    ph = lon(:);

    % Preallocate result matrices
    num_coeffs = (max_deg + 1)^2 - 1;
    A_r = zeros(num_pts, num_coeffs);
    A_th = zeros(num_pts, num_coeffs);
    A_ph = zeros(num_pts, num_coeffs);

    sin_th = sin(th);
    cos_th = cos(th);
    
    idx = 0;
    for n = 1:max_deg
        % Set radial scaling based on source type
        if strcmpi(source_type, 'ext')
            r_scale = r.^(n-1);
            radial_factor = -n;
        else
            r_scale = r.^(-(n+2));
            radial_factor = (n+1);
        end

        % Compute Schmidt-normalized Legendre polynomials
        Pnm_all = legendre(n, cos_th, 'sch')';
        
        % Compute derivatives (dPnm) using recurrence relations
        dPnm = zeros(num_pts, n+1);
        dPnm(:, 1) = -sqrt(n*(n+1)/2) * Pnm_all(:, 2); % m=0
        dPnm(:, n+1) = sqrt(n/2) * Pnm_all(:, n);      % m=n
        
        if n > 1
            dPnm(:, 2) = (sqrt(2*n*(n+1)) * Pnm_all(:, 1) - sqrt((n+2)*(n-1)) * Pnm_all(:, 3)) / 2;
        elseif n == 1
            dPnm(:, 2) = sqrt(2) * dPnm(:, 2);
        end
        
        for m = 2:n-1
            dPnm(:, m+1) = (sqrt((n+m)*(n-m+1)) * Pnm_all(:, m) - sqrt((n+m+1)*(n-m)) * Pnm_all(:, m+2)) / 2;
        end

        % Fill matrices for each order m
        for m = 0:n
            c_lon = cos(m*ph);
            s_lon = sin(m*ph);
            Pnm = Pnm_all(:, m+1);
            D_Pnm = dPnm(:, m+1);

            if m == 0
                idx = idx + 1;
                A_r(:, idx)  = radial_factor .* r_scale .* Pnm;
                A_th(:, idx) = -r_scale .* D_Pnm;
                % A_ph remains 0 for m=0
            else
                % g_n^m (Internal) or q_n^m (External)
                idx = idx + 1;
                A_r(:, idx)  = radial_factor .* r_scale .* c_lon .* Pnm;
                A_th(:, idx) = -r_scale .* c_lon .* D_Pnm;
                A_ph(:, idx) = (m * r_scale .* s_lon .* Pnm) ./ sin_th;

                % h_n^m (Internal) or s_n^m (External)
                idx = idx + 1;
                A_r(:, idx)  = radial_factor .* r_scale .* s_lon .* Pnm;
                A_th(:, idx) = -r_scale .* s_lon .* dPnm(:, m+1);
                A_ph(:, idx) = (-m * r_scale .* c_lon .* Pnm) ./ sin_th;
            end
        end
    end
end