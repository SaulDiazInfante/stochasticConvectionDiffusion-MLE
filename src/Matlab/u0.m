function u0xy = u0(x, y, Lx, Ly)
%U0 Initial scalar field used as initial condition
%   u0xy = u0(x,y,Lx,Ly) returns the value of the initial condition at
%   coordinates x and y. Inputs x and y may be scalars or arrays of the
%   same size. The function is vectorized and uses elementwise operations.

    % Validate inputs
    if ~(isscalar(x) || isscalar(y) || isequal(size(x), size(y)))
        error('u0:InvalidInput', 'x and y must be scalars or arrays of the same size.');
    end

    frac = 1/4;

    % Normalized radius (elementwise operations)
    xr = (x - Lx/4) ./ (frac * Lx);
    yr = (y - 0.6 * Ly) ./ (frac * Ly);
    rad = sqrt(xr .^ 2 + yr .^ 2);

    % Initialize output and apply compact support formula where rad <= 1
    u0xy = zeros(size(rad));
    mask = (rad <= 1);
    if any(mask(:))
        u0xy(mask) = (1 + cos(pi * rad(mask))) / 2;
    end
end



