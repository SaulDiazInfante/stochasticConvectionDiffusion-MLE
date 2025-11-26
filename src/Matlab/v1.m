function v1xy = v1(x, y, L1, L2)
    % First component of the velocity field
    %
    % Computes the first component of the velocity field at point (x, y).
    % The velocity field is a polynomial term modulated by a Gaussian envelope.
    %
    % Inputs:
    %   x: x-coordinate (scalar or array)
    %   y: y-coordinate (scalar or array)
    %   L1: domain length in x-direction
    %   L2: domain length in y-direction
    %
    % Output:
    %   v1xy: velocity field value at (x, y)
    
    % Domain parameters
    r = 1;
    
    % Gaussian envelope centered at (L1/2, L2/2)
    gaussian_envelope = exp(-0.6 * ((x - L1/2)^2 + (y - L2/2)^2)^r);
    
    % Polynomial velocity component
    polynomial_component = y + x + 2*cos(x/2) + sin(x/2);
    
    % Combined velocity field
    v1xy = polynomial_component * gaussian_envelope;
end