function v2xy = v2(x, y, L1, L2)
    % Second component of the velocity field
    %
    % Computes the second component of the velocity field at point (x, y).
    % The velocity field is a polynomial term modulated by a Gaussian envelope.
    %
    % Inputs:
    %   x: x-coordinate (scalar or array)
    %   y: y-coordinate (scalar or array)
    %   L1: domain length in x-direction
    %   L2: domain length in y-direction
    %
    % Output:
    %   v2xy: velocity field value at (x, y)
    
    % Envelope parameters
    r = 1;
    
    % Gaussian envelope centered at (L1/2, L2/2)
    gaussian_envelope = exp(-0.6 * ((x - L1/2)^2 + (y - L2/2)^2)^r);
    
    % Polynomial velocity component
    polynomial_component = ...
        (-8/17) * exp(x/2) * cos(2*x) + (2/17) * exp(x/2) * sin(2*x) - y;
    
    % Combined velocity field
    v2xy = polynomial_component * gaussian_envelope;
end