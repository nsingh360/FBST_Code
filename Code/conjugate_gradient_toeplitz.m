function [x, iter, rsold] = conjugate_gradient_toeplitz(T, b, tol, max_iter)
    % Conjugate Gradient Method for Hermitian Semi-Positive Definite Toeplitz System
    %
    % T: Toeplitz matrix (Hermitian semi-positive definite)
    % b: right-hand side vector
    % tol: tolerance for convergence
    % max_iter: maximum number of iterations
    
    % Ensure the matrix is Hermitian
    assert(isequal(T, T'), 'Matrix T must be Hermitian.');
    
    n = length(b);      % Dimension of the system
    x = zeros(n, 1);    % Initial guess
    r = b - toep_mult(T, x);  % Initial residual
    p = r;              % Initial search direction
    rsold = r' * r;     % Initial residual norm squared
    
    for iter = 1:max_iter
        Ap = toep_mult(T, p);  % Efficient multiplication with Toeplitz matrix
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        
        if sqrt(rsnew) < tol
            break;
        end
        
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end
