function X = proj_col21(V, z)

    norms_V = sqrt(sum(V.^2));
    n = length(norms_V);

    lambda_max = max(norms_V);
    lambda_min = -z/n;


    % Find the zero of f(\lambda) using bisection
    MAX_ITER = 25;
    for i = 1:MAX_ITER
        lambda = (lambda_max + lambda_min) / 2;
        if auxiliary(norms_V, lambda, z) > 0
            lambda_min = lambda;
        else
            lambda_max = lambda;
        end
    end

    alpha = (1 - (lambda ./ norms_V));
    alpha(alpha < 0) = 0;

    X = V * diag(alpha);


function f = auxiliary(norms_V, lambda, z)

    f = sum(norms_V((norms_V - lambda) >= 0) - lambda) - z;

