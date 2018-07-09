load_data = 1;
save_data = 0;
 
ensemble = 'laplace';

delta_list = [.02 .1:.1:.9];
n_list = [200];

if load_data
    fname = sprintf('Example_nonGaussian_%s.mat', ensemble);
    load(fname);
else
    n_iter = 100;
    frob2_pinvA = zeros(length(n_list), length(delta_list), n_iter);
    frob2_pinvA_nrm = zeros(length(n_list), length(delta_list), n_iter);
    frob2_spinvA = zeros(length(n_list), length(delta_list), n_iter);
    frob2_spinvA_nrm = zeros(length(n_list), length(delta_list), n_iter);

    i = 0;
    for n = n_list
        i = i + 1;
        j = 0;
        for delta = delta_list
            j = j + 1;
            for iter = 1:n_iter
                fprintf('n %d, delta %f, iter %d\n', n, delta, iter);
                m = round(delta*n + 1);
                
                switch ensemble
                    % normalize all to unit variance
                    case 'student'
                        A = trnd(3, m, n) / sqrt(3);
                    case 'gaussian'
                        A = randn(m, n);
                    case 'bernoulli'
                        A = sign(randn(m, n));
                    case 'uniform'
                        A = sqrt(12) * (rand(m, n) - 0.5);
                    case 'laplace'
                        A = laprnd(m, n, 0, 1);
                    otherwise
                        error('Unknown matrix ensemble')
                end
                
                pinvA = pinv(A);
                spinvA = ginvadmm(A, 'l1', 0);
                frob2_pinvA(i, j, iter) = norm(pinvA, 'fro')^2;
                frob2_pinvA_nrm(i, j, iter) = n/m*norm(pinvA, 'fro')^2;
                frob2_spinvA(i, j, iter) = norm(spinvA, 'fro')^2;
                frob2_spinvA_nrm(i, j, iter) = n/m*norm(spinvA, 'fro')^2;
            end
        end
    end

    if save_data
        fname = sprintf('Example_nonGaussian_%s.mat', ensemble);
        save(fname);
    end
end

%% curves

figure(1); clf;
hold all;
delta_cont = .02:.01:.9;
astar2lim_pinv = 1 ./ (1 - delta_cont);
t_1 = sqrt(2)*erfcinv(delta_cont);
astar2lim_spinv = 1 ./ (sqrt(2/pi).*exp(-t_1.^2/2).*t_1 - delta_cont.*t_1.^2) - 1./delta_cont;

plot(delta_cont, astar2lim_spinv, 'linewidth', 4);
scatter(delta_list, mean(frob2_spinvA_nrm(1, :, :), 3), 120, 's', 'linewidth', 2);

% figure(2); clf;
% hold all;
plot(delta_cont, astar2lim_pinv, 'linewidth', 4);
scatter(delta_list, mean(frob2_pinvA_nrm(1, :, :), 3), 120, 'd', 'linewidth', 2);

xlabel('Delta');
ylabel('Scaled squared Frobenius norm');

set(gca, 'FontSize', 20);

axis tight



