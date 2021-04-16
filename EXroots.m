%%
%1
a = 1;
f = @(x)x.*x-a;
x = fzero(f,0.5)
%%
%2a
clear
a = -3; b = 3;
m = 100; x = linspace(a,b,m);
y = x.^2+2*x-1-sin(x);
plot (x, y);
grid on
xlabel('x'); ylabel('y');
z = ginput(1);
[k,l] = fzero(@(x)x.^2+2*x-1-sin(x), z(1));
hold on
plot(k,l,'r*',z(1),z(2),'g*')
hold off
%%
%2b
syms x
Eq = x*x+2*x-1-sin(x);
a = solve(Eq)
%%
%3a
clear
x = -3:.01:3;
y = sin(exp(x));
plot(x,y);
grid on
z = ginput(1);
[k,l] = fzero(@(x)sin(exp(x)), z(1));
hold on
plot(k,l,'r*',z(1),z(2),'g*')
hold off
%%
%3b
clear
x = 0:.01:2*pi;
y = sin(x.*(1-x));
plot(x,y);
grid on
z = ginput(1);
[r,l] = fzero(@(x)sin(x.*(1-x)), z(1));
hold on
plot(r,l,'r*',z(1),z(2),'g*')
hold off
%%
%3c
clear
x = 0:.01:4*pi;
y = x.*sin(x)-cos(x);
plot(x,y);
grid on
z = ginput(1);
[k,l] = fzero(@(x)x.*sin(x)-cos(x), z(1));
hold on
plot(k,l,'r*',z(1),z(2),'g*')
hold off
%% 
%3d
clear
x = pi/2:.01:3*pi;
y = (sin(x)).^2+(0.5-1./x).*cos(x)-0.5;
plot(x,y);
grid on
z = ginput(1);
[m,l] = fzero(@(x)(sin(x)).^2+(0.5-1./x).*cos(x)-0.5, z(1));
hold on
plot(m,l,'r*',z(1),z(2),'g*')
hold off
%% 
%3e
clear
x = -2*pi:.01:6*pi;
y = 5.*exp(-0.1*x).*sin(x)-0.1.*x;
plot(x,y);
grid on
z = ginput(1);
[n,l] = fzero(@(x)5.*exp(-0.1*x).*sin(x)-0.1.*x, z(1));
hold on
plot(n,l,'r*',z(1),z(2),'g*')
hold off
%%
%3f
clear
n = input("Введите степень Х в уравнении x^n-1=0");
x = -1.5:.01:1.5;
y = x.^n-1;
plot(x,y);
grid on
z = ginput(1);
[p,l] = fzero(@(x)x.^n-1, z(1));
hold on
plot(p,l,'r*',z(1),z(2),'g*')
hold off
%%
%4
syms u 
e = u^2 + 1;
S1 = solve(e == 0, u);
disp(S1);
S2 = fzero(@(v)v^2 + 1, 1);
disp(S2)
%%
%5
x = linspace(-5,5,100);
y = cos(x) - exp(0.001+x.*x)
plot(x,y)
syms u 
e = cos(u) - exp(0.001 + u^2) ;
S1 = solve(e, u);
disp(S1);
S2 = fzero(@(v)cos(v) - exp(0,001+v.^2), 1);
disp(S2)
%% Задача 1а[Метод Ньютона]
clear
a = -10; b = 10; eps = 0.001; h = 0.1; iter = 100;
x = linspace(a,b,100);
y = 1 + x.*sin(x);
plot(x,y);
hold on;
plot(x, y*0);
z = ginput(1);
z = z(1);
for i = 1:iter
    line([z,z],[0, 1+z*sin(z)], 'Color', 'blue');
    x = z-1:0.1:z+1;
    y = (x-z)*((z+h)*sin(z+h)-(z)*sin(z))/h+1+z*sin(z);
    plot(x,y,"--" ,'Color', 'red');
    k = z - h*(1+z*sin(z))/((z+h)*sin(z+h)-(z)*sin(z));
    plot(z,1+z*sin(z),'.');
    plot(k, 1+k*sin(k), '.');
    z = k;
    if abs(1+k*sin(k)) < eps
        break;
    end
end
disp(k);
hold off
%% 
%Задача 1а[Метод деления отрезка пополам]
clear
a = -10; b = 10; eps = 0.001;
x = linspace(a,b,100);
y = 1 + x.*sin(x);
plot(x,y); 
hold on;
plot(x, y*0);
grid on;
[z1, z2] = ginput(2);
z1 = z1(1); z2 = z2(1);
if (sign(1 + z1.*sin(z1)) == sign(1 + z2.*sin(z2)) && sign(1 + ((z1+z2)/2).*sin((z1+z2)/2)) == sign(1 + z2.*sin(z2)) && sign(1 + z1.*sin(z1)) == sign(1 + ((z1+z2)/2).*sin((z1+z2)/2)))
    plot(z1,1 + z1.*sin(z1),'.');
    plot((z1+z2)/2,1 + ((z1+z2)/2).*sin((z1+z2)/2),'.');
    plot(z2,1 + z2.*sin(z2),'.');
    disp("На этом отрезке нет корней")
else
    for i = 1:iter
        plot(z1,1 + z1.*sin(z1),'.');
        plot((z1+z2)/2,1 + ((z1+z2)/2).*sin((z1+z2)/2),'.');
        plot(z2,1 + z2.*sin(z2),'.');
        x = (z1+z2)/2;
        if abs(1 + x.*sin(x)) < eps
            break;
        end
        if sign(1 + x.*sin(x)) == sign(1 + z2.*sin(z2))
            z2 = x;
        else
            z1 = x;
        end
    end
    disp(x);
    disp(1 + x.*sin(x));
    disp(i);
end
hold off