% Definimos la cantidad de puntos
aux = 30;
% Hacemos nuestra partición
X = linspace(0,1,aux+2);
ceros = zeros(1,aux+2);
Ti0 = zeros(1,aux+2);
Ti1 = zeros(1,aux+2);
meh = zeros(1,aux+2);
% Definimos nuestras funciones
for i = 1:aux+2
    Ti0(i) = sin(2*pi*X(i));
    Ti1(i) = X(i)*(1-X(i));
end
% Graficamos nuestros resultados.
close all;
T = Temperatura(aux,aux,Ti1,Ti0,ceros,ceros);
[XX, YY] = meshgrid(X,X);
surf(XX,YY,T);
colormap(jet(5));
colorbar;
xlabel('Eje X','fontweight','bold','Fontsize',10);
ylabel('Eje Y','fontweight','bold', 'Fontsize',10);
zlabel('Temperatura','fontweight','bold', 'Fontsize',10);
title('Distribución de calor en placa metálica', 'Fontsize',12);