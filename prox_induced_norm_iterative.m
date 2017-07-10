function [Xhat] = prox_induced_norm_iterative(V,lambda)
%
% This function computes the proximal operator for the l1 induced norm of
% a matrix that is, the solution to the following optimization problem:
%
% minimize_X |X|_1 + 1/2/lambda * trace[(X-V)'(X-V)]
%
% where |X|_1 = max_i |xi|_1, X = [x1,...,xm] and |xi| = sum_j abs(xi(j)).
%

% get data dimensions
[n,m] = size(V);

% variables
U = abs(V);
u = sum(U);

% compute initial t
t = (sum(u) - n*lambda)/m;
% set of columns involved
M = find( t < u );
Mold = M;
N = length(M);

% search loop
while true
  
  % recompute J
  J = true(n,N);
  Jold = J;
  
  % inner search loop
  while true
    
    % compute t
    t = 0; w = 0;
    for ii = 1:N
      t = t + mean(U(J(:,ii),M(ii)));
      w = w + 1/nnz(J(:,ii));
    end
    t = ( t - lambda ) / w ;
    
    % compute mu
    mu = zeros(1,N);
    for ii = 1:N
      mu(ii) = ( sum( U( J(:,ii), M(ii) )) - t )/nnz(J(:,ii));
    end
    
    % compute
    J = bsxfun(@gt,U(:,M),mu);
    
    if all(Jold==J)
      break;
    end
    
    % update set
    Jold = J;
    
  end
  
  % update M
  M = find( t < u );
  N = length(M);
  
  % exit condition
  if numel(M) == numel(Mold)
    if Mold == M
      break;
    end
  end
  
  % store past value of M
  Mold = M;
  
end

% output
Xhat = U;
Xhat(:,M) = max(bsxfun(@minus,Xhat(:,M),mu),0);
negindx = V < 0;
Xhat(negindx) = -Xhat(negindx);
