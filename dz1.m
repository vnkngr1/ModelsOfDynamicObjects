% --- Основной скрипт моделирования ---

function dYdt = ode_fun(t, Y, k)
    % Параметры системы для Варианта 4
    m = 20;
    l0 = 1; % Недеформированная длина пружины, м
    g = 9.81; % Ускорение свободного падения, м/с^2

    % Внешняя сила g(t)
    g_t = 9.815 * sin(0.05 * pi * t);
    
    % Обобщенные координаты
    theta = Y(1); % y1
    theta_dot = Y(2); % y2
    Delta = Y(3); % y3
    Delta_dot = Y(4); % y4
    
    % Система ДУ (4 уравнения первого порядка)
    dYdt(1,1) = theta_dot; % d(theta)/dt = theta_dot (y2)
    
    % d(theta_dot)/dt = y2_dot
    % theta_dot_dot = - (2*Delta_dot*theta_dot)/(l0 + Delta) - ((g + g_t)/(l0 + Delta)) * sin(theta)
    theta_dot_dot = - (2 * Delta_dot * theta_dot) / (l0 + Delta) - ((g + g_t) / (l0 + Delta)) * sin(theta);
    dYdt(2,1) = theta_dot_dot; 
    
    dYdt(3,1) = Delta_dot; % d(Delta)/dt = Delta_dot (y4)
    
    % d(Delta_dot)/dt = y4_dot
    % Delta_dot_dot = (l0 + Delta)*theta_dot^2 - (k/m)*Delta + (g + g_t)*cos(theta)
    Delta_dot_dot = (l0 + Delta) * theta_dot^2 - (k/m) * Delta + (g + g_t) * cos(theta);
    dYdt(4,1) = Delta_dot_dot;
end

% Параметры
k_nom = 5; % Номинальный k
k_max = 6.0; % +20%
k_min = 4.0; % -20%
k_values = [k_min, k_nom, k_max];

% Начальные условия (произвольные)
Y0 = [0.5; 0; 0.5; 0]; % [theta0; theta_dot0; Delta0; Delta_dot0]

% Вектор времени
t_max = 100; % c
dt = 0.01;   % c
tspan = 0:dt:t_max;

results = cell(1, 3);
titles = {'k = 4.0 (-20%)', 'k = 5.0 (Номинал)', 'k = 6.0 (+20%)'};

for i = 1:length(k_values)
    k = k_values(i);
    % Решение системы ДУ
    [T, Y] = ode45(@(t, Y) ode_fun(t, Y, k), tspan, Y0);
    results{i}.T = T;
    results{i}.Y = Y;
end

% --- Построение графиков ---

% Графики изменения обобщенных координат во времени
figure(1);
sgtitle('Влияние k на Динамику Системы (k = 4.0, 5.0, 6.0)', 'FontSize', 14);

% График для theta(t)
subplot(2, 1, 1);
hold on;
for i = 1:length(results)
    plot(results{i}.T, results{i}.Y(:, 1), 'DisplayName', titles{i});
end
hold off;
title('Изменение Угла \theta(t)');
xlabel('Время t, с');
ylabel('Угол \theta, рад');
legend('Location', 'best');
grid on;

% График для Delta(t)
subplot(2, 1, 2);
hold on;
for i = 1:length(results)
    plot(results{i}.T, results{i}.Y(:, 3), 'DisplayName', titles{i});
end
hold off;
title('Изменение Удлинения Пружины \Delta(t)');
xlabel('Время t, с');
ylabel('Удлинение \Delta, м');
legend('Location', 'best');
grid on;

% Фазовые портреты
figure(2);
sgtitle('Фазовые Портреты Системы (k = 4.0, 5.0, 6.0)', 'FontSize', 14);

% Фазовый портрет для \theta (theta - theta_dot)
subplot(1, 2, 1);
hold on;
for i = 1:length(results)
    plot(results{i}.Y(:, 1), results{i}.Y(:, 2), 'DisplayName', titles{i});
end
hold off;
title('Фазовый Портрет: \theta - \dot{\theta}');
xlabel('\theta, рад');
ylabel('\dot{\theta}, рад/с');
legend('Location', 'best');
grid on;

% Фазовый портрет для \Delta (Delta - Delta_dot)
subplot(1, 2, 2);
hold on;
for i = 1:length(results)
    plot(results{i}.Y(:, 3), results{i}.Y(:, 4), 'DisplayName', titles{i});
end
hold off;
title('Фазовый Портрет: \Delta - \dot{\Delta}');
xlabel('\Delta, м');
ylabel('\dot{\Delta}, м/с');
legend('Location', 'best');
grid on;