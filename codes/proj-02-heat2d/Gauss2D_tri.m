function [xi, eta, weight] = Gauss2D_tri(n_int)

if n_int == 3
    % 3-point rule
    xi = [1/6, 2/3, 1/6];
    eta = [1/6, 1/6, 2/3];
    weight = [1/6, 1/6, 1/6] * 0.5;

else
    error('Unsupported order for Gauss quadrature in triangle.');
end
end
