%% Zadanie_1
%Aproksymacja sinusa oraz cosiunusa
%Potrzebne do rachunków w metodzie epsilonów
clear;
syms x;

subplot(2,1,1);
fplot(cos(x),[0,1],'r');
hold on
fplot(1 - x^2/2,[0,1],'b');
hold off
grid on;
xlabel('x')

subplot(2,1,2);
fplot(sin(x),[0,1],'r');
hold on
fplot(x,[0,1],'b');
hold off
grid on;
xlabel('x')
%% Zadanie_1 błędy 
%Wykresy zależności błędów od x
clear;
syms x;

T(x) = 2*x^2 * cot(x^2) + (x/((x+2)*log(x+2))) - 3*x^3;
K1(x) = x^2 * cot(x^2);
K2(x) = 1/log(x+2);
K3(x) = -x^3;
K4 = 1;
K6 = -1;

subplot(3,2,1);
fplot(T,[0,1],'r');
xlabel('x');
ylabel("T(x)");
grid on;

subplot(3,2,2);
fplot(K1,[0,1],'r');
xlabel('x');
ylabel("K1(x)");
grid on;

subplot(3,2,3);
fplot(K2,[0,1],'r');
xlabel('x');
ylabel("K2(x)");
grid on;

subplot(3,2,4);
fplot(K3,[0,1],'r');
xlabel('x');
ylabel("K3(x)");
grid on;

subplot(3,2,5);
fplot(K4,[0,1],'r');
xlabel('x');
ylabel("K4(x),K5(x),K7(x),K8(x)");
grid on;

subplot(3,2,6);
fplot(K6,[0,1],'r');
xlabel('x');
ylabel("K6(x)");
grid on;

%% Zadanie_2 metoda maksymalizacji sumy modułów
clear;
N = 1000;
x = linspace(0,1,N);

T = 2.*(x.^2).*cot(x.^2)+(x./((x+2).*log(x+2)))-3.*x.^3;
K1 = (x.^2).*cot(x.^2);
K2 = 1./log(x+2);
K3 = -x.^3;
K4 = 1;
K5 = 1;
K6 = -1;
K7 = 1;
K8 = 1;
eps1 = 4*10^(-13);
delta_err = abs(T) + abs(K1) + abs(K2) + abs(K3) + abs(K4) + abs(K5) + abs(K6) + abs(K7) + abs(K8);

%Obliczanie maksymalnej sumy modułów
max_delta_err = max(delta_err);

%błąd całkowity wyznaczania wartości y
sup_delta_err_zadanie_2 = max_delta_err*eps1

%% Zadanie_3 metoda symulacyjna
clear;
N = 1000;
x = linspace(0,1,N);
eps1 = 4*10^(-13);

%Określenie ile będzie wszystkich możliwych błędów 
%(9 pozycji na każdej eps lub -eps czyli 2^9)
dimension = 2^9;
M = ([1,dimension]);

y = (sin(x.^2).*log(x+2))./exp(x.^3);
for i=1:dimension
    k = int2bit(i,9);
    k(k==1)=eps1;
    k(k==0)=(-1).*eps1;
    
    y_err = (sin((x.*(1+k(1))).^2*(1+k(2))).*log((x.*(1+k(1))+2).*(1+k(3))).*(1+k(5)).*(1+k(6)).*(1+k(8)).*(1+k(9)))./(exp(((x.*(1+k(1))).^3).*(1+k(4))).*(1+k(7)));
    delta_err = abs((y_err-y)./y);
    M(1,i) = max(delta_err);
end
sup_delta_err_zadanie_3 = max(M)
