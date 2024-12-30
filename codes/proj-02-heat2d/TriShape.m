function [val, grad_xi_eta] = TriShape(aa, xi, eta)
    if aa == 1
        val = 1 - xi - eta;
        grad_xi_eta = [-1; -1];
    elseif aa == 2
        val = xi;
        grad_xi_eta = [1; 0];
    elseif aa == 3
        val = eta;
        grad_xi_eta = [0; 1];
    else
        error('Node index should be 1, 2, or 3');
    end
end