load('result_smoother');

subplot(2, 2, 1);
errorbar(mus_f(1, :), sigs_f(1, :)); hold on;
errorbar(mus_fb(1, :), sigs_fb(1, :)); hold on;
xlabel('Time');
ylabel('x');
legend('Filter', 'Smoother');
hold off;

subplot(2, 2, 2);
errorbar(mus_f(2, :), sigs_f(2, :)); hold on;
errorbar(mus_fb(2, :), sigs_fb(2, :)); hold on;
xlabel('Time');
ylabel('y');
legend('Filter', 'Smoother');
hold off;

subplot(2, 2, 3);
errorbar(mus_f(3, :), sigs_f(3, :)); hold on;
errorbar(mus_fb(3, :), sigs_fb(3, :)); hold on;
xlabel('Time');
ylabel('vx');
legend('Filter', 'Smoother');
hold off;

subplot(2, 2, 4);
errorbar(mus_f(4, :), sigs_f(4, :)); hold on;
errorbar(mus_fb(4, :), sigs_fb(4, :)); hold on;
xlabel('Time');
ylabel('vy');
legend('Filter', 'Smoother');
hold off;

