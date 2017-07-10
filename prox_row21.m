function X = prox_row21(V, lambda)

normr_V = sqrt(sum(abs(V).^2, 2));
scale = (1 - lambda ./ normr_V);
scale(scale < 0) = 0;

X = diag(scale) * V;
