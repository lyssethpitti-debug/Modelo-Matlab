%% =========================================================
%  ECUACION DE TRANSPORTE 1D POR DIFERENCIAS FINITAS
%  u_t + c u_x = 0
%
%  Este script permite:
%  1) Comparar derivadas numericas:
%     - diferencia hacia adelante
%     - diferencia hacia atras
%     - diferencia centrada
%
%  2) Resolver la ecuacion de transporte usando esquema upwind
%
%  3) Elegir distintas condiciones iniciales:
%     - trivial
%     - senoidal
%     - escalon
%     - gaussiana
%
%  Autor: [Tu nombre]
%  Curso: Metodos Numericos
%% =========================================================

clear; clc; close all;

%% =========================================================
%  PARTE 0. PARAMETROS GENERALES
%% =========================================================

% Velocidad de transporte
c = 1;

% Dominio espacial: 0 <= x <= 1
a = 0;
b = 1;

% Tiempo final
T = 0.5;

% Paso espacial y temporal
h = 0.05;
k = 0.02;

% Numero de Courant lambda = c*k/h
lambda = c*k/h;

fprintf('=============================================\n');
fprintf('PARAMETROS DEL PROBLEMA\n');
fprintf('c       = %.4f\n', c);
fprintf('h       = %.4f\n', h);
fprintf('k       = %.4f\n', k);
fprintf('lambda  = %.4f\n', lambda);
fprintf('=============================================\n\n');

% Verificacion de estabilidad del esquema upwind
if lambda < 0 || lambda > 1
    warning('La condicion CFL no se cumple: 0 <= lambda <= 1');
    warning('El metodo puede ser inestable o poco confiable.');
else
    fprintf('La condicion CFL se cumple: 0 <= lambda <= 1\n\n');
end

%% =========================================================
%  PARTE 1. MALLA ESPACIAL Y TEMPORAL
%% =========================================================

% Vector de nodos espaciales
x = a:h:b;

% Vector de tiempos
t = 0:k:T;

% Numero de nodos espaciales y temporales
Nx = length(x);
Nt = length(t);

fprintf('Numero de nodos espaciales Nx = %d\n', Nx);
fprintf('Numero de nodos temporales Nt = %d\n\n', Nt);

%% =========================================================
%  PARTE 2. COMPARACION DE DERIVADAS NUMERICAS
%
%  Usamos una funcion de prueba conocida:
%      f(x) = sin(pi*x)
%
%  Su derivada exacta es:
%      f''(x)? NO
%      f'(x) = pi*cos(pi*x)
%
%  Aqui compararemos:
%      adelante: (f(i+1)-f(i))/h
%      atras:    (f(i)-f(i-1))/h
%      centrada: (f(i+1)-f(i-1))/(2h)
%% =========================================================

fprintf('=============================================\n');
fprintf('PARTE 2: COMPARACION DE DERIVADAS NUMERICAS\n');
fprintf('Funcion de prueba: f(x) = sin(pi*x)\n');
fprintf('Derivada exacta:   f''(x) = pi*cos(pi*x)\n');
fprintf('=============================================\n\n');

% Definimos la funcion de prueba
f  = @(x) sin(pi*x);
df = @(x) pi*cos(pi*x);

% Evaluamos funcion y derivada exacta en la malla
fx        = f(x);
df_exacta = df(x);

% Reservamos memoria para aproximaciones numericas
dif_adelante = NaN(1, Nx);
dif_atras    = NaN(1, Nx);
dif_centrada = NaN(1, Nx);

% Reservamos memoria para errores absolutos
err_adelante = NaN(1, Nx);
err_atras    = NaN(1, Nx);
err_centrada = NaN(1, Nx);

% ---------------------------------------------------------
% Diferencia hacia adelante
% Se puede calcular desde i = 1 hasta Nx-1
% Porque necesitamos f(x_{i+1})
% ---------------------------------------------------------
for i = 1:Nx-1
    dif_adelante(i) = (fx(i+1) - fx(i))/h;
    err_adelante(i) = abs(df_exacta(i) - dif_adelante(i));
end

% ---------------------------------------------------------
% Diferencia hacia atras
% Se puede calcular desde i = 2 hasta Nx
% Porque necesitamos f(x_{i-1})
% ---------------------------------------------------------
for i = 2:Nx
    dif_atras(i) = (fx(i) - fx(i-1))/h;
    err_atras(i) = abs(df_exacta(i) - dif_atras(i));
end

% ---------------------------------------------------------
% Diferencia centrada
% Se puede calcular desde i = 2 hasta Nx-1
% Porque necesitamos f(x_{i-1}) y f(x_{i+1})
% ---------------------------------------------------------
for i = 2:Nx-1
    dif_centrada(i) = (fx(i+1) - fx(i-1))/(2*h);
    err_centrada(i) = abs(df_exacta(i) - dif_centrada(i));
end

% ---------------------------------------------------------
% Construccion de tabla
% ---------------------------------------------------------
TablaDerivadas = table( ...
    x', fx', df_exacta', ...
    dif_adelante', dif_atras', dif_centrada', ...
    err_adelante', err_atras', err_centrada', ...
    'VariableNames', { ...
    'x', 'f_x', 'derivada_exacta', ...
    'dif_adelante', 'dif_atras', 'dif_centrada', ...
    'error_adelante', 'error_atras', 'error_centrada'});

disp(TablaDerivadas);

% Errores maximos ignorando NaN
max_err_adelante = max(err_adelante, [], 'omitnan');
max_err_atras    = max(err_atras, [], 'omitnan');
max_err_centrada = max(err_centrada, [], 'omitnan');

fprintf('\nErrores maximos de las aproximaciones:\n');
fprintf('Error maximo adelante = %.6f\n', max_err_adelante);
fprintf('Error maximo atras    = %.6f\n', max_err_atras);
fprintf('Error maximo centrada = %.6f\n\n', max_err_centrada);

%% =========================================================
%  PARTE 3. ELECCION DE CONDICION INICIAL
%
%  Opciones:
%   1 -> Trivial:   u(x,0)=0
%   2 -> Senoidal:  u(x,0)=sin(pi*x)
%   3 -> Escalon:   u(x,0)=1 en [0.2,0.6], 0 en otro caso
%   4 -> Gaussiana: u(x,0)=exp(-alpha*(x-x0)^2)
%
%  La condicion de frontera se fija como:
%      u(0,t)=0
%
%  Esto es coherente con c > 0, pues la informacion entra por x=0.
%% =========================================================

fprintf('=============================================\n');
fprintf('PARTE 3: SELECCION DE CONDICION INICIAL\n');
fprintf('1 -> Trivial\n');
fprintf('2 -> Senoidal\n');
fprintf('3 -> Escalon\n');
fprintf('4 -> Gaussiana\n');
fprintf('=============================================\n');

opcion = input('Seleccione una opcion (1,2,3,4): ');

% Vector de condicion inicial
u0 = zeros(1, Nx);

switch opcion
    
    case 1
        % -------------------------------------------------
        % Condicion inicial trivial
        % u(x,0) = 0
        % -------------------------------------------------
        u0 = zeros(1, Nx);
        nombre_ci = 'Trivial';
        
    case 2
        % -------------------------------------------------
        % Condicion inicial senoidal
        % u(x,0) = sin(pi*x)
        % -------------------------------------------------
        u0 = sin(pi*x);
        nombre_ci = 'Senoidal';
        
    case 3
        % -------------------------------------------------
        % Condicion inicial tipo escalon
        % u(x,0) = 1 si 0.2 <= x <= 0.6
        %          0 en otro caso
        % -------------------------------------------------
        for i = 1:Nx
            if x(i) >= 0.2 && x(i) <= 0.6
                u0(i) = 1;
            else
                u0(i) = 0;
            end
        end
        nombre_ci = 'Escalon';
        
    case 4
        % -------------------------------------------------
        % Condicion inicial gaussiana
        % u(x,0)=exp(-alpha*(x-x0)^2)
        % -------------------------------------------------
        alpha = 100;
        x0 = 0.4;
        u0 = exp(-alpha*(x - x0).^2);
        nombre_ci = 'Gaussiana';
        
    otherwise
        error('Opcion no valida. Debe elegir 1, 2, 3 o 4.');
end

fprintf('\nCondicion inicial seleccionada: %s\n\n', nombre_ci);

%% =========================================================
%  PARTE 4. RESOLUCION DE LA ECUACION DE TRANSPORTE
%
%  Ecuacion:
%      u_t + c u_x = 0
%
%  Esquema upwind explicito para c > 0:
%
%      u_i^{n+1} = (1-lambda)u_i^n + lambda*u_{i-1}^n
%
%  donde:
%      lambda = c*k/h
%
%  Frontera:
%      u(0,t)=0
%% =========================================================

fprintf('=============================================\n');
fprintf('PARTE 4: RESOLUCION DE LA ECUACION DE TRANSPORTE\n');
fprintf('Esquema upwind explicito\n');
fprintf('=============================================\n\n');

% Matriz de solucion:
% filas    -> tiempo
% columnas -> espacio
U = zeros(Nt, Nx);

% ---------------------------------------------------------
% Imponemos la condicion inicial
% Primera fila corresponde a t = 0
% ---------------------------------------------------------
U(1, :) = u0;

% ---------------------------------------------------------
% Imponemos la condicion de frontera:
% u(0,t)=0 para todo t
% La primera columna corresponde a x = 0
% ---------------------------------------------------------
U(:, 1) = 0;

% ---------------------------------------------------------
% Aplicamos el esquema upwind
% n: indice temporal
% i: indice espacial
%
% IMPORTANTE:
% Para cada punto interior usamos:
% U(n+1, i) = (1-lambda)*U(n,i) + lambda*U(n,i-1)
% ---------------------------------------------------------
for n = 1:Nt-1
    for i = 2:Nx
        U(n+1, i) = (1 - lambda)*U(n, i) + lambda*U(n, i-1);
    end
end

fprintf('Calculo de la solucion completado.\n\n');

%% =========================================================
%  PARTE 5. TABLA PARCIAL DE RESULTADOS DE LA SOLUCION
%
%  Mostraremos una tabla con algunos tiempos:
%      t=0
%      tiempo intermedio
%      tiempo final
%% =========================================================

indice_medio = round(Nt/2);

TablaSolucion = table( ...
    x', U(1,:)', U(indice_medio,:)', U(end,:)', ...
    'VariableNames', {'x', 'u_t0', 'u_tmedio', 'u_tfinal'});

fprintf('=============================================\n');
fprintf('TABLA PARCIAL DE LA SOLUCION NUMERICA\n');
fprintf('=============================================\n\n');

disp(TablaSolucion);

%% =========================================================
%  PARTE 6. GRAFICAS
%% =========================================================

% ---------------------------------------------------------
% Grafica 1: Comparacion de derivadas
% ---------------------------------------------------------
figure;
plot(x, df_exacta, 'k-', 'LineWidth', 2); hold on;
plot(x, dif_adelante, 'o-', 'LineWidth', 1.2);
plot(x, dif_atras, 's-', 'LineWidth', 1.2);
plot(x, dif_centrada, 'd-', 'LineWidth', 1.2);
grid on;
xlabel('x');
ylabel('Derivada');
title('Comparacion de derivadas numericas');
legend('Exacta', 'Hacia adelante', 'Hacia atras', 'Centrada', ...
       'Location', 'best');

% ---------------------------------------------------------
% Grafica 2: Errores de las derivadas
% ---------------------------------------------------------
figure;
plot(x, err_adelante, 'o-', 'LineWidth', 1.2); hold on;
plot(x, err_atras, 's-', 'LineWidth', 1.2);
plot(x, err_centrada, 'd-', 'LineWidth', 1.2);
grid on;
xlabel('x');
ylabel('Error absoluto');
title('Errores de las aproximaciones de la derivada');
legend('Error adelante', 'Error atras', 'Error centrada', ...
       'Location', 'best');

% ---------------------------------------------------------
% Grafica 3: Condicion inicial
% ---------------------------------------------------------
figure;
plot(x, U(1,:), 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('u(x,0)');
title(['Condicion inicial seleccionada: ', nombre_ci]);

% ---------------------------------------------------------
% Grafica 4: Evolucion temporal de la solucion
% ---------------------------------------------------------
figure;
plot(x, U(1,:), 'LineWidth', 2); hold on;
plot(x, U(indice_medio,:), 'LineWidth', 2);
plot(x, U(end,:), 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('u(x,t)');
title(['Evolucion de la solucion - CI: ', nombre_ci]);
legend(['t = ', num2str(t(1))], ...
       ['t = ', num2str(t(indice_medio))], ...
       ['t = ', num2str(t(end))], ...
       'Location', 'best');

% ---------------------------------------------------------
% Grafica 5: Mapa de calor de la solucion
% ---------------------------------------------------------
figure;
imagesc(x, t, U);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('x');
ylabel('t');
title(['Mapa de calor de la solucion - CI: ', nombre_ci]);

%% =========================================================
%  PARTE 7. MENSAJES INTERPRETATIVOS
%% =========================================================

fprintf('=============================================\n');
fprintf('INTERPRETACION DE RESULTADOS\n');
fprintf('=============================================\n');

fprintf('\n1) En la tabla de derivadas numericas:\n');
fprintf('   - La diferencia centrada suele ser mas precisa.\n');
fprintf('   - Las diferencias hacia adelante y hacia atras\n');
fprintf('     suelen tener mayor error.\n');

fprintf('\n2) En la ecuacion de transporte:\n');
fprintf('   - El esquema upwind es apropiado para c > 0.\n');
fprintf('   - La informacion se desplaza hacia la derecha.\n');
fprintf('   - La frontera u(0,t)=0 hace que no entre nueva informacion.\n');

fprintf('\n3) Sobre la condicion inicial elegida:\n');
switch opcion
    case 1
        fprintf('   - La solucion permanece nula en todo el dominio.\n');
    case 2
        fprintf('   - La senal senoidal se transporta hacia la derecha.\n');
    case 3
        fprintf('   - El escalon permite observar difusion numerica.\n');
    case 4
        fprintf('   - La gaussiana permite observar transporte de un pulso suave.\n');
end

fprintf('\n4) Si deseas experimentar:\n');
fprintf('   - Cambia h y k para ver el efecto de la malla.\n');
fprintf('   - Cambia c para ver la velocidad de transporte.\n');
fprintf('   - Cambia la condicion inicial para comparar perfiles.\n');
fprintf('   - Verifica siempre que 0 <= lambda <= 1.\n\n');