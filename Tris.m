function [x] = Tris(U,b)
% Solución del sistema triangular superior: U * x = b, por medio de 
% sustitución hacia atrás. La matriz U es triangular superior

% 'n' representa la cantidad de variables de nuestro sistema.
n = length(b);
% x es el vector donde guardaremos la resolución del sistema.
x = zeros(1,n);
% El ciclo for nos ayuda a realizar la sustitución hacia atras,
% despejando una variable en cada iteración.
for i = n : -1 : 1
    suma = 0;
    for z = i+1 : n
        suma = suma + U(i,z)*x(z);
    end
    x(i) = (b(i)-suma)/U(i,i);
end
end