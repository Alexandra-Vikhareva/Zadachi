%% задание 1
P = [1 -3.55 5.1 -3.1]
x = linspace(-10,10,1000);
F = polyval(P,x);
plot(x, [F;0*F]);
grid on;
a = roots(P)
k=0;
for i=1:length(a)
    if (isreal(a(i))==1)
        k=k+1;
        b(k)=a(i)
    end
end
hold on
for i=1:length(k)
    plot(b(i),0,"r*")
end
hold off
%% задание 2 А
P = [1 0.1 0.2 -0.2 -2 1]
a = roots(P)
k=0; r = 0;
for i=1:length(a)
    if (isreal(a(i))==1)
        k=k+1;
        b(k)=a(i)
    end
    if ((a(i))^2>r^2)
        r = a(i)
    end
end
m = max(b)
n = min(b)
x = linspace(n,m,100);
F = polyval(P,x);
hold on
plot(x,[F; 0*F]);
K=[5 0.4 0.6 -0.4 -2];
l = m/2;
k1=polyval(K,l);
k01=polyval(P,l);
Y1 = k01+k1.*(x-l);
plot(x,Y1);
l1=-m/2
k2=-1/polyval(K,l1);
k02=polyval(P,l1);
l2=l1+cos(k2);
k002=k02+sin(k2);
line([l1,l2],[k02,k002])
th = 0:0/001:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);
plot(xunit, yunit);
hold off
%% задание 3
u1 = [2 -3 4 -5 6]; v1 = [1 -3 1];
u2 = [1 -3 -1 -1]; v2 = [3 -2 1];
[q1,r1]=deconv(u1, v1);
[q2,r2]=deconv(u2,v2);
disp(q1, r1, q2, r2)
%% задание 4 А
clc;
clear;
P = [1 -2 6 -10 16];
x0 = 4;
P1(1) = P(1);
for i = 2:5
    P1(i) = P1(i-1)*x0 + P(i)
end

if (P1(5) == polyval(P, 4))
    disp('значение полинома в точке x = 4 :'); disp(P1(5));
end;

for i = 1:4
    P2(i) = P1(i)
end
P3(1) = P2(1);
for i = 2:4
    P3(i) = P3(i-1)*x0 + P2(i)
end

if (P3(4) == polyval(polyder(P), 4))
    disp('значение производной полинома в точке x = 4 :'); disp(P3(4));
end;
%% задание 4 Б
clc;
clear;
P = [1 1+2*i 0 -(1+3*i) 0 7];
x0 = -2-i;
P1(1) = P(1);
for k = 2:6
    P1(k) = P1(k-1)*x0 + P(k)
end

if (P1(6) == polyval(P, -2-i))
    disp('значение полинома в точке x = -2-i :'); disp(P1(6));
end;

for k = 1:5
    P2(k) = P1(k)
end
P3(1) = P2(1);
for k = 2:5
    P3(k) = P3(k-1)*x0 + P2(k)
end

if (P3(5) == polyval(polyder(P), -2-i))
    disp('значение производной полинома в точке x = -2-i :'); disp(P3(5));
end;

%% задание 5 А
clear
e=1*10^(-5);
n = input("Введите степень N в уравнении X^2*N-N*X^(N+1)+N*X^(N-1)-1");
P(2*n+1)=-1;
P(1)=1;
for i=2*n+1:-1:1
    if i == n+2
        P(i)=P(i)+n
    end
    if i == n
        P(i)=P(i)-n
    end
end
a = roots(P)
for i=1:(2*n)
    if (abs(imag(a(1))))<e
        d(i)=real(a(i))
    end
end
m=100;
l=-100;
for i=1:length(a)
    if m>a(i)
        m=a(i)
    end
    if l<a(i)
        l=a(i)
    end
end
x = linspace(m-1,l+1,100);
y = polyval(P,x);
plot(x,y)
k = 1;
b(1)=d(1);
for i=1:(2*n)
    u=0;
    for j =1:k
        if imag(d(i))==imag(b(j)) && real(d(i))==real(b(j))
            u=u+1
        end
    end
    if (u==0)
        k=k+1;
        b(k)=d(i)
    end
end
r=0;
for i=1:k
    u=0;
    for j=1:(2*n)
        if imag(b(i))==imag(d(j)) && real(b(i))==real(d(j))
            u=u+1
        end
    end
    if u>1
        r=r+1
        c(r)=b(i)
    end
end
hold on
for i=1:r
    plot(c,0,'*r')
end
hold off
%% задание 6 А
B = [1 0 0];
A = [1 4 1 -6];
[r3,p3,K3]=residue(B,A)
%% задание 6 Б
B = [1 3];
A = [1 -1 1 -1];
[r4,p4,K4]=residue(B,A)
%% задание 6 В
B = [1 0 0];
A = [1 0 0 0 -1];
[r5,p5,K5]=residue(B,A)
%% задание 7 А
x = linespace(-5,5,100);
y = 1/.x;
plot(x,y)
disp("Корней нет")
%% задание 7 Б
A = [1 -1 1 -1 1; 0 0 0 0 1; 1 1 1 1 1; 16 8 4 2 1; 81 27 9 3 1]
B = [6; 5; 0; 3; 2]
x = inv(A)*B;
P = x';
a = roots(P);
m = max(a);
n = min(a);
X = linspace(n,m,100);
Y = polyval(P,X)
plot(X,Y)
%% задание 8 А
P = [1 0 -7 -12 6 36];
r = roots(P)
l = length(r);
m = 0;
for i = 1:l
    if abs(P(i))>=m
        m = abs(P(i))
    end
end
n = 1 + abs(m/P(1));
x = linspace(-n, n, 100);
y = x.^5 - 7.*x.^3 -12.*x.^2-26.*x+12;
plot(x,y)
grid on;
disp(r)
%% задание 8 Б
P = [1 -6 15 -14];
r = roots(P)
l = length(r);
m = 0;
for i = 1:l
    if abs(P(i))>=m
        m = abs(P(i))
    end
end
n = 1 + abs(m/P(1));
x = linspace(-n, n, 100);
y = x.^3-6.*x.^2+15.*x-14;
plot(x,y)
grid on
disp(r)
%% задание 8 В
P = [1 -5 11 -1 -18 20 -8];
r = roots(P)
l = length(r);
m = 0;
for i = 1:l
    if abs(P(i))>=m
        m = abs(P(i))
    end
end
n = 1 + abs(m/P(1));
x = linspace(-n, n, 100);
y = x.^6 -6.*x.^5+11.*x.^4-x.^3-18.*x.^2+20.*x-8;
plot(x,y)
grid on
disp(r)

%% задание 9
n = randi(10)
A = randi(n,n)
e = eig(A)
P = poly(e)
disp(P)
if (poly(A)==P)
    disp ("Матрица является корнем полинома")
else 
    disp("Матрица НЕ является корнем полинома")
end
%% задание 10 А
clear
syms x
P = x.^4+2.*x.^3-x.^2-4.*x-2;
Q = x.^4+x.^3-x.^2-2.*x-2;
[G,S,T] = gcd(P,Q)
disp('НОД:'); disp(G)
disp('Линейное представление: '); disp(S); disp(T)
%% задание 10 Б
clear
syms x
P = x.^5+3.*x.^4+x.^3+x.^2+3.*x+1;
Q = x.^4+2.*x.^3+x+2;
[G,S,T] = gcd(P,Q)
disp('НОД:'); disp(G)
disp('Линейное представление: '); disp(S); disp(T)
%% задание 10 В
clear
syms x
P = 4.*x.^4-2.*x.^3-16.*x.^2+5.*x+9;
Q = 2.*x.^3-x.^2-5.*x+4;
[G,S,T] = gcd(P,Q)
disp('НОД:'); disp(G)
disp('Линейное представление: '); disp(S); disp(T)
%% задание 11
A = ones(100,1)
x=linspace(0,2,101)
P = poly(A)
plot(x,P)