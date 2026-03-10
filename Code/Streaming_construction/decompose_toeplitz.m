function [T1, T2, error] = decompose_toeplitz(A,D,max_iter, tol)
    % Input:
    % A - matrix to decompose
    % max_iter - maximum number of iterations
    % tol - tolerance for convergence
    
    % Output:
    % T1, T2 - Toeplitz matrices such that A ≈ T1 * T2
    % error - Frobenius norm of the approximation error
    
    [m, n] = size(A);
    
    % Initialize T1 and T2 as random Toeplitz matrices
    c1 = D(:,1); 
    r1 = D(1,:);
    T1 = toeplitz(c1, r1);

    Er = inv(D)*A;
    c2 = Er(:,1); 
    r2 = Er(1,:);
    T2 = toeplitz(c2, r2);

    % Iterative updates
    for iter = 1:max_iter
        
        % Update T2 given T1
        T2 = T1 \ A;
        % Ensure T2 is Toeplitz
        T2 = toeplitz(T2(:,1), T2(1,:));
        
        % Calculate the approximation error
        current_error = norm(A - T1 * T2, 'fro');
        
        % Display the iteration number and current error
        fprintf('Iteration %d: Error = %f\n', iter, current_error);
        
        % Check for convergence
        if current_error < tol
            break;
        end
    end

    % Return the final error
    error = current_error;
end
