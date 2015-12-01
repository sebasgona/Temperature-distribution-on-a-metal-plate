function [T] = Temperatura(n,m,T0j,T1j,Ti0,Ti1)
% Se calcula la temperatura en la placa metálica con una partición
% en el eje X de n + 2 puntos igualmente espaciados y con m + 2 puntos
% sobre el eje Y.
% T0j es un vector de dimensión m + 2 donde T0j(k) = T0, k
% T1j es un vector de dimensión m + 2 donde T1j(k) = T1, k
% T i0 es un vector de dimensión n + 2 donde T i0(k) = Tk, 02
% T i1 es un vector de dimensión n + 2 donde T i1(k) = Tk, 1
% La salida, T es una matriz de orden (n + 2)x(m + 2) tal que T(i, j) es
% la temperatura en el nodo (xi,yj), 0 ? i ? n + 1, 0 ? j ? m + 1.
T = zeros(n+2,m+2);
AUX = zeros(n*m);
temp = zeros(1,(n*m));

% Asignamos los valores de temperatura ya conocidos a los lados
T(1:n+2,1) = Ti0(1:n+2);
T(1:n+2,n+2) = Ti1(1:n+2);
T(1,1:m+2) = T0j(1:m+2);
T(n+2,1:m+2) = T1j(1:m+2);
% Ponemos especial enfasis en las esquinas ya que contienen dos fuentes.
T(1,1) = Ti0(1)+T0j(1);
T(1,n+2) = Ti0(1)+T1j(n+2);
T(n+2,1) = Ti1(n+2)+T0j(1);
T(n+2,n+2) = Ti1(n+2)+T1j(n+2);

% Utilizamos la variable ec para manejar el renglon donde insertaremos la
% ecuación de temperatura del nodo (i,j).
ec = 1;
cant = n+1;
for i=2 : cant
    for j=2 : cant
% Ponemos siempre el coeficiente 4 en la diagonal ya que ese es el nodo que
% estamos analizando para obtener su temperatura. 
% 4T(i,j) = T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1)
        AUX(ec,ec) = 4;
        
% Checamos cada caso posible y llenamos los coeficientes de la ecuación
% respectiva. (Esquinas, lados, puntos interiores.)
        if i==2 || j==2
            if i==2 && j==2
                AUX(ec,ec+1) = -1;
                AUX(ec,ec+n) = -1;
                temp(ec) = T(1,2)+T(2,1);
                elseif i==2 && j==cant
                        AUX(ec,ec-1) = -1;
                        AUX(ec,ec+n) = -1;
                        temp(ec) = T(1,j)+T(i,cant+1);
                    elseif i==cant && j==2
                        AUX(ec,ec+1) = -1;
                        AUX(ec,ec-n) = -1;
                        temp(ec) = T(cant,1)+T(cant+1,2);
                        elseif i == 2
                            AUX(ec,ec+1) = -1;
                            AUX(ec,ec-1) = -1;
                            AUX(ec,ec+n) = -1;
                            temp(ec) = T(1,j);
            else
                AUX(ec,ec+1) = -1;
                AUX(ec,ec-n) = -1;
                AUX(ec,ec+n) = -1;
                temp(ec) = T(i,1);
            end
        elseif i==cant || j==cant
                if i==cant && j==cant
                    AUX(ec,ec-1) = -1;
                    AUX(ec,ec-n) = -1;
                    temp(ec) = T(cant+1,cant)+T(cant,cant+1);
                    elseif i==cant
                        AUX(ec,ec+1) = -1;
                        AUX(ec,ec-n) = -1;
                        AUX(ec,ec-1) = -1;
                        temp(ec) = T(cant+1,j);
                else
                    AUX(ec,ec-1) = -1;
                    AUX(ec,ec-n) = -1;
                    AUX(ec,ec+n) = -1;
                    temp(ec) = T(i,cant+1);
                end
        else
            if ec -1 > 0
            AUX(ec,ec-1) = -1;
            end
            if ec+1 < (n*m)^2
            AUX(ec,ec+1) = -1;
            end
            if ec+n <= (n*m)^2
            AUX(ec,ec+n) = -1;
            end
            if ec-n <= n*m
            AUX(ec,ec-n) = -1;
            end
            temp(ec) = 0;
        end
        ec = ec+1;
    end
end

% Obtenemos nuestra factorización LU y nuestra matríz pivote
[P,L,U] = MiLuPiv(AUX);
% Pivoteamos los resultados del vector temp para que este asociado a la
% factorización LU
temp = P'*temp';
% Resolvemos el sistemas utilizando nuestra factorización.
y = Trin(L,temp);
temp = Tris(U,y);

% Llenamos nuestra matríz T con los datos obtenidos.
pos = 1;
for i=2 : cant
    for j=2 : cant
        T(i,j) = temp(pos);
        pos = pos+1;
    end
end
end