function X = prox_l1l2(V, lambda)

X = V - lambda*proj_col21(V, lambda);
