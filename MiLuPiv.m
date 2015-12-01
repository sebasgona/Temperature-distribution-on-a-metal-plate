function [P,L,U] = MiLuPiv(A)
% Factorizacón LU con pivoteo parcial de la matriz cuadrada A.
% Matriz de permutaciones en P. Matriz triangular inferior con unos en la
% diagonal en L. Matriz triangular superior en U. Se tiene que P*A = L*U.

[n,n] = size(A);
P = eye(n);
L = eye(n);
AUX = zeros(1,n);
U = zeros(n);

for i=1 : n-1
% Revisamos que la posición en la diagonal sea distinta de cero, en caso 
% contrario buscamos el valor más grande en la columna e intercambiamos los
% renglones.
    if A(i,i) == 0
        MAX = A(i,i+1);
        pos = i+1;
        for k=(i+2) : n
            if A(k,i) > MAX
                MAX = A(k,i);
                pos = k;
            end
        end
        AUX = A(pos,1:n);
        A(pos,1:n) = A(i,1:n);
        A(i,1:n) = AUX;
% Guardamos la permutación realizada dentro de P. 
        AUX = P(pos,1:n);
        P(pos,1:n) = P(i,1:n);
        P(i,1:n) = AUX;
    end
% Volvemos cero los elementos debajo del elemento i-ésimo sobre la dagonal.
    L(i+1 : n ,i) = A(i+1 : n , i)/A(i,i);
    A(i+1:n , i+1 : n) = A(i+1:n , i+1 : n) -  L(i+1 : n,i)*A(i,i+1:n);
    A(i+1:n,i) = zeros(n-i,1);
end
U = A;
end
