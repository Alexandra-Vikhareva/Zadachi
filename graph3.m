%% пример 2(в полярных координатах)
f = 0:.01:2*pi;
r = sin(f)./f;
polar(f,r)

%% пример 2(в дпск)
x = @(fi)cos(fi).*sin(fi)/fi;
y = @(fi)(sin(fi))^2/fi;
ezplot(x,y,[0, 2*pi]);

%% пример 3
x = 0:.01:2;
y = humps(x);
plot(x, y)
[a,b] = max(y);
b = (b-1)/100;
hold on
plot(b, a , "Color", "r", "Marker","o");
c = find(y<40 & y>20);
c = (c-1)./100;
d = humps(c);
plot(c, d, ".", "Marker", "o", "Color", "r");
hold off

%% задание 1
x = -1:.01:-0.01;
x2 = 0.01:.01:1;
y = (1+x.^3)./(x.^2);
y2 =  (x2.^3 +1)./(x2.^2);
plot(x,y, x2, y2, "Color", "b")
title("график X+1/X^2")

%% задание 7
x = -5:.01:5;
y = x.*sin(pi*x);
plot(x,y)
title("график X*sin(pi*x)")
axis

%% задание 8
a = 100;
b = 2;
n = 4;
k = 14;
fi = -pi/2:.01:3*pi/2;
r = a./(a+(fi-pi/2)).^n.*(b-sin(k.*fi)-cos(m.*fi))
polar(fi,r)

%% задание 12
a = 1/4;
b = 1/16;
m = 8;
n = 8;
k = 1;
x = @(t)cos(t)-a*cos(m.*t)+b*sin(n.*t);
y = @(t)sin(t)-a*sin(m.*t)+b*cos(n.*t);
ezplot(x, y, [-5, 5]);

%% задание 13.1
x = [-3:.01:-0.01 0.01:.01:0.99 1.01:.01:1.99 2.01:.01:5];
y = 1./x+1./(x-1)+1./(x-2);
plot(x,y)

%% задание 13.2
a=7;
b=3;
fi=-pi:.01:pi;
rho = a*cos(fi)-b;
polar(fi, rho)

%% задание 13.3
n = 2;
m = 2;
fi = -5:.01:5;
x = cos(n.*fi).*(cos(fi)).^m;
y = sin(n.*fi).*(sin(fi)).^m;
plot(x,y)

%% задание 14.1 цепная линия
a = 4;
x = -5:.01:5
y = (exp(x./a)+exp(-x./a))*a/2;
plot(x,y)

%% задание 14.2 спираль Ферма
a = 2;
fi = -5*pi:.01:5*pi;
rho1 = sqrt(a^2*fi);
rho2=-sqrt(a^2*fi);
hold on
polar(fi, rho1) 
polar(fi, rho2)
hold off

%% задание 14.3 лист Хабенихта
a = 8;
theta=-pi:.01:pi;
rho=1+7*cos(a.*theta)+4*(sin(a.*theta)).^2+3*(sin(a.*theta)).^4;
polar(theta, rho)
