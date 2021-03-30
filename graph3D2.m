%% Пример 2
t = 0:.1:6*pi;
x = sin(t);
y = 2*cos(t);
plot3(x,y,t, "b")
line ([-1,-1],[-21.9911485751,18.0088514249],[0,20])
%% Пример 
[X,Y] = meshgrid(-pi:.1:pi, -pi:.1:pi);
Z = 20 -X.^2 -Y.^2;
surf(X,Y,Z);
line ([0,0],[-1,-3],[19,20]);
R = 0.1;
fi = 0:.1:6*pi;
r = R*fi;
x = r.*cos(fi);
y = r.*sin(fi);
z = 20 -x.^2 -y.^2;
hold on;
plot3(x,y,z);
z = 0.*z;
plot3(x,y,z);
hold off
%% Задание 1
a = 1;
b = 2;
c = 3;
d = 4;
[X,Y]=meshgrid(1:.5:10,1:.5:10);
Z = -(d/c) - (b/c).*Y - (a/c).*X;
surf(X,Y,Z)
%% Задание 2.A
a = 5;
[X,Y] = meshgrid(-5:.1:5, -5:.1:5);
Z = (a .* sin ((X.^2+Y.^2).^(1/2)))./((X.^2+Y.^2).^(1/2));
subplot(2,1,1)
mesh(X,Y,Z);
subplot(2,1,2)
plot3(X,Y,Z)
%% Задание 2.Б
[X,Y] = meshgrid(-5:.1:5, -5:.1:5);
Z = -X.*sin(X)-Y.*cos(Y);
subplot(2,1,1)
mesh(X,Y,Z);
subplot(2,1,2)
plot3(X,Y,Z)
%% Задание 3.1(a)
a = 5;
[X,Y] = meshgrid(-5:0.25:5,-5:0.25:5);
Z = a.*X.*exp(-X.^2-Y.^2);
subplot(3,1,1);
surf(X,Y,Z);
subplot(3,1,2);
surfc(X,Y,Z);
subplot(3,1,3);
mesh(X,Y,Z);
%% Задание 3.2(b)
a = 2;
c = 5;
alpha1 = 0; alpha2 = 5; beta1 = 0; beta2 = 2*pi;
N1 = 40; N2 = 40;
alpha = linspace(alpha1, alpha2, N1);
beta = linspace(beta1, beta2, N2);
[Alpha, Beta] = meshgrid(alpha, beta);
x =a .* cosh(Alpha).*cos(Beta);
y =a .* cosh(Alpha).*sin(Beta);
z =c .* sinh(Alpha);
figure
subplot(3,1,1);
surf(x,y,z);
subplot(3,1,2);
surfc(x,y,z);
subplot(3,1,3);
mesh(x,y,z);
%% Задание 3.3(d)
a = 2;
c = 5;
alpha1 = 0; alpha2 = 5; beta1 = 0; beta2 = 2*pi;
N1 = 40; N2 = 40;
alpha = linspace(alpha1, alpha2, N1);
beta = linspace(beta1, beta2, N2);
[Alpha, Beta] = meshgrid(alpha, beta);
x =a .* sinh(Alpha).*cos(Beta);
y =a .* sinh(Alpha).*sin(Beta);
z =c .* cosh(Alpha);
figure
subplot(3,1,1);
surf(x,y,z);
hold on
surf(x,y,-z);
hold off
subplot(3,1,2);
surfc(x,y,z);
hold on
surfc(x,y,-z);
hold off
subplot(3,1,3);
mesh(x,y,z);
hold on
mesh(x,y,-z);
hold off