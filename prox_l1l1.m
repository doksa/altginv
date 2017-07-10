function X = prox_l1l1(V, lambda)

X = ( V' - lambda*projL1Inf(V'/lambda, 1, ones(size(V', 1))) )';
