function [x] = Trin(L,b)
% Soluci�n del sistema triangular inferior: L * x = b por medio de 
% sustituci�n hacia adelante. La matriz L es triangular inferior.

% 'n' representa la cantidad de variables de nuestro sistema.
n = length(b);
% x es el vector donde guardaremos la resoluci�n del sistema.
x = zeros(1,n);
% El ciclo for nos ayuda a realizar la sustituci�n hacia adelante,
% despejando una variable en cada iteraci�n.
for i = 1 : n
    suma = 0;
    for z = 1 : i-1
        suma = suma + L(i,z)*x(z);
    end
    x(i) = (b(i)-suma)/L(i,i);
end
end

