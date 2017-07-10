function X = ginvadmm(A, norm_name, proj)

% X = ginvadmm(A, norm_name, proj)
%
% Computes a generalized inverse of A minimizing norm_name using the 
% alternating direction method of multipliers (ADMM).
%
% INPUT:  A         ... input matrix to be inverted
%
%         norm_name ... 'l1'   -> entrywise l1 norm
%                       'l1l1' -> induced l1->l1 norm
%                       'l1l2' -> induced l1->l2 norm
%
%         proj      ... 0      -> minimize the norm of X
%                       1      -> minimize the norm of X*A
%
% OUTPUT: X         ... the sought generalized inverse
%
% By Ivan Dokmanic & Remi Gribonval 2014


% Set up some termination criteria
MAX_ITER = 1000;
RELTOL   = 1e-5;

[m, n] = size(A);

% These things probably only work for fat matrices
if m > n
    error('A should be fat!');
end

% Select the correct proximal operator
if strcmp(norm_name, 'l1')
    proxfcn = @prox_l1;
elseif strcmp(norm_name, 'l1l1')
    proxfcn = @prox_l1l1;
elseif strcmp(norm_name, 'l1l2')
    proxfcn = @prox_l1l2;
elseif strcmp(norm_name, 'row21')
    proxfcn = @prox_row21;
end


N      = null(A);
Proj_N = N*N';
B      = pinv(A);

if proj
    U = zeros(n);
    Z = zeros(n);
    X = zeros(n, m);

    lambda = 0.5;
    if strcmp(norm_name, 'row21')
        mu = lambda / norm(A, 'fro') / (n/2); % Dirty, but needed for U([0, 1]) input matrices with n>>m!
    else
        mu = lambda / norm(A, 'fro') / (n/4); % You need to make tables of these values...
    end

    for k = 1:MAX_ITER
        % Update X, Z, and U
        X_old = X;
        X = B + Proj_N*(X - mu/lambda*(X*A - Z + U)*A');
        Z_old = Z;
        Z = feval(proxfcn, X*A + U, lambda);
        U = U + X*A - Z;
        
        % Fix this
        if norm(Z_old - Z, 'fro') / norm(Z_old, 'fro') < RELTOL
            break;
        end
    end
else
    U = zeros(n, m);
    Z = zeros(n, m);
    X = zeros(n, m);

    lambda = .1;
    for k = 1:MAX_ITER
        % Update X, Z, and U
        X = B + Proj_N*(Z - U);
        Z_old = Z;
        Z = feval(proxfcn, X + U, lambda);

        % Fix this
        if norm(Z_old - Z, 'fro') / norm(Z_old, 'fro') < RELTOL
            break;
        end

        U = U + X - Z;
    end
end

if k == MAX_ITER
    fprintf('Warning: MAXITER reached. Relative error is %f. Norm is %s.\n', norm(Z_old - Z, 'fro') / norm(Z_old, 'fro'), norm_name);
end


