function plot_results(type, result, pos_vec, hist_times)

switch type

    case 'sis'

        figure;
        plot(result.tau1, result.tau2, 'b-', 'LineWidth', 1.5); hold on;
        plot(pos_vec(1,:), pos_vec(2,:), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8);
        xlabel('x_1');
        ylabel('x_2');
        title('SIS estimated trajectory');
        legend('estimated trajectory','base stations','Location','best');
        grid on;
        axis equal;

        figure;
        plot(0:result.m, result.ESS, 'LineWidth', 1.5);
        xlabel('n');
        ylabel('ESS');
        title('ESS');
        grid on;

        for j = 1:length(hist_times)
            figure;
            histogram(result.histW(:,j), 50);
            xlabel('normalized weight');
            ylabel('frequency');
            title(['Weights at n = ', num2str(hist_times(j))]);
            grid on;
        end

    case 'sisr'

        figure;
        plot(result.tau1, result.tau2, 'b-', 'LineWidth', 1.5); hold on;
        plot(pos_vec(1,:), pos_vec(2,:), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8);
        xlabel('x_1');
        ylabel('x_2');
        title('SISR estimated trajectory');
        legend('estimated trajectory','base stations','Location','best');
        grid on;
        axis equal;

        figure;
        stairs(0:result.m, result.cmd_mode, 'LineWidth', 1.2);
        yticks(1:5);
        yticklabels({'stay','east','north','south','west'});
        xlabel('n');
        ylabel('most probable command');
        title('Most probable driving command');
        grid on;

        figure;
        plot(0:result.m, result.cmd_prob, 'LineWidth', 1.2);
        xlabel('n');
        ylabel('posterior probability');
        title('Posterior probabilities of commands');
        legend('stay','east','north','south','west','Location','best');
        grid on;

    case 'calibration'

        figure;
        plot(result.varsigma_grid, result.loglik_grid, 'o-', 'LineWidth', 1.5);
        xlabel('\varsigma');
        ylabel('normalized log-likelihood estimate');
        title('Approximate normalized log-likelihood');
        grid on;

        figure;
        plot(result.tau1_hat, result.tau2_hat, 'b-', 'LineWidth', 1.5); hold on;
        plot(pos_vec(1,:), pos_vec(2,:), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8);
        xlabel('x_1');
        ylabel('x_2');
        title(['Estimated trajectory under \varsigma hat = ', ...
            num2str(result.varsigma_hat)]);
        legend('estimated trajectory','base stations','Location','best');
        grid on;
        axis equal;

end

end