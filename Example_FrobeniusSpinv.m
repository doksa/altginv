load_data = 0;

delta_list = .1:.1:.9;
n_list = [50, 100, 200, 500, 1000];

if load_data
    load('Example_FrobeniusSpinv.mat');
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
                A = randn(m, n);
                pinvA = pinv(A);
                spinvA = ginvadmm(A, 'l1', 0);
                frob2_pinvA(i, j, iter) = norm(pinvA, 'fro')^2;
                frob2_pinvA_nrm(i, j, iter) = n/m*norm(pinvA, 'fro')^2;
                frob2_spinvA(i, j, iter) = norm(spinvA, 'fro')^2;
                frob2_spinvA_nrm(i, j, iter) = n/m*norm(spinvA, 'fro')^2;
            end
        end
    end

    save Example_FrobeniusSpinv
end


%%
delta_cont = delta_list;
astar2lim_pinv = 1 ./ (1 - delta_cont);
t_1 = sqrt(2)*erfcinv(delta_cont);
astar2lim_spinv = 1 ./ (sqrt(2/pi).*exp(-t_1.^2/2).*t_1 - delta_cont.*t_1.^2) - 1./delta_cont;

d_sel = [1, 2, 4, 7];
n_sel = 2:5;

figure(1); clf;
hold all;
for i = n_sel
    for j = d_sel
        re = frob2_pinvA_nrm(i, j, :);
        re = re(:); 
        scatter(n_list(i)*ones(size(re)), re, 'filled');
        scatter(n_list(i), mean(re), 100, 'ks', 'linewidth', 2);
    end
end

for j = d_sel
    plot([0, 1000], [astar2lim_pinv(j), astar2lim_pinv(j)], 'k-.', 'linewidth', 1);
end

xlabel('n');
ylabel('Frobenius norm');

figure(2); clf;
hold all;
for i = n_sel
    for j = d_sel
        re = frob2_spinvA_nrm(i, j, :);
        re = re(:); 
        scatter(n_list(i)*ones(size(re)), re, 'filled');
        scatter(n_list(i), mean(re), 100, 'ks', 'linewidth', 2);
    end
end

for j = d_sel
    plot([0, 1000], [astar2lim_spinv(j), astar2lim_spinv(j)], 'k-.', 'linewidth', 1);
end

xlabel('n', 'fontsize', 14);
ylabel('Frobenius norm', 'fontsize', 14);

%% boxplots

figure(3); clf;
hold all;

k = 0;
for j = d_sel
    delta_list(j)
    plot([0, 1000], [astar2lim_spinv(j), astar2lim_pinv(j)], 'k-.', 'linewidth', 1);

    mat_j = reshape(frob2_spinvA_nrm(:, j, :), length(n_list), n_iter);
    boxplot(mat_j', 'widths', .5);
    
    ylim([2, 6]);
    k = k + 1;
end


%% curves

figure(3); clf;
hold all;
delta_cont = .05:.01:.9;
astar2lim_pinv = 1 ./ (1 - delta_cont);
t_1 = sqrt(2)*erfcinv(delta_cont);
astar2lim_spinv = 1 ./ (sqrt(2/pi).*exp(-t_1.^2/2).*t_1 - delta_cont.*t_1.^2) - 1./delta_cont;

plot(delta_cont, astar2lim_spinv);
scatter(delta_list, mean(frob2_spinvA_nrm(4, :, :), 3), 120, 's', 'linewidth', 1);
scatter(delta_list, mean(frob2_spinvA_nrm(3, :, :), 3), 120, 's', 'linewidth', 1);
scatter(delta_list, mean(frob2_spinvA_nrm(2, :, :), 3), 120, 's', 'linewidth', 1);
scatter(delta_list, mean(frob2_spinvA_nrm(5, :, :), 3), 200, 'sk', 'linewidth', 2);

plot(delta_cont, astar2lim_pinv);
scatter(delta_list, mean(frob2_pinvA_nrm(4, :, :), 3), 120, 'd', 'linewidth', 1);
scatter(delta_list, mean(frob2_pinvA_nrm(3, :, :), 3), 120, 'd', 'linewidth', 1);
scatter(delta_list, mean(frob2_pinvA_nrm(2, :, :), 3), 120, 'd', 'linewidth', 1);
scatter(delta_list, mean(frob2_pinvA_nrm(5, :, :), 3), 200, 'ok', 'linewidth', 2);

xlabel('Delta', 'fontsize', 14);
ylabel('Scaled squared Frobenius norm', 'fontsize', 14);


