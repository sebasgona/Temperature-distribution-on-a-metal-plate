function [x] = Trin(L,b)
% Solución del sistema triangular inferior: L * x = b por medio de 
% sustitución hacia adelante. La matriz L es triangular inferior.

% 'n' representa la cantidad de variables de nuestro sistema.
n = length(b);
% x es el vector donde guardaremos la resolución del sistema.
x = zeros(1,n);
% El ciclo for nos ayuda a realizar la sustitución hacia adelante,
% despejando una variable en cada iteración.
for i = 1 : n
    suma = 0;
    for z = 1 : i-1
        suma = suma + L(i,z)*x(z);
    end
    x(i) = (b(i)-suma)/L(i,i);
end
end

