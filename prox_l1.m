function X = prox_l1(V, lambda)

X = wthresh(V, 's', lambda);