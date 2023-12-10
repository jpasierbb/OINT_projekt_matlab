%% Zadanie 1 - 
clear;
R = 1000;
x = linspace(-1,1, R);
f = (2/3).*cos(pi.*(x-1/3)).*exp(x);

N = 30; %Zmieniając N zmieniamy liczbę próbek do aproksymacji
f_x = linspace(-1,1,N);
f_k = (2/3).*cos(pi.*(f_x-1/3)).*exp(f_x);

figure(1);
plot(x ,f, 'r');
hold on;
plot(f_x, f_k, 'bo');
hold off;
title(sprintf('Zadanie 1, N = %d',N));
xlabel('x');
ylabel('f(x)');
grid on;
legendstr = sprintf('Próbki do aproksymacji dla N = %d',N);
legend('Orygnalna funkcja f(x)', legendstr, "Location", "northwest")

%% Zadanie 2
% Dane f_x i f_k
clear;
x = linspace(-1,1, 1000);
K = 50; %Ilosc funkcji do aproksymacji, musi byc mnijesza od ilosci wezlow
N = 100; %Ilosc wezlow
x_k = linspace(-1,1,K);
%[M,Y,P] = rysowanie(N,K);
%plot(linspace(-1,1, 20),M)
%Q_1000 = aprox(1000,K);
[Q_1000, Q_N] = aprox(N,K);
figure(2);
zad1(N);
hold on;
plot(x,Q_1000, "b");
plot(linspace(-1,1,N),Q_N, "ok")
hold off;
title(sprintf('Zadanie 2, K = %d oraz N = %d', K, N));
xlabel('x');
ylabel('f(x)');
grid on;
legend('Orygnalna funkcja f(x)', 'Aproksymacja', "Location", "northwest")

%% Zadanie 3
clear;
err = zeros(46,46);
err_max = zeros(46,46);
for i=5:50
    for j=4:(i-1)
        err(i-4,j-3) = error(i,j);
        err_max(i-4,j-3) = error_max(i,j);
    end
end
figure(3);
mesh(linspace(5,50,46), linspace(4,49,46), log10(err));
title("Zadanie 3 - Wskaźniki błędu średniokwadratowego");
xlabel('N');
ylabel('K');
zlabel('err');

figure(4);
mesh(linspace(5,50,46), linspace(4,49,46), log10(err_max));
title("Zadanie 3 - Wskaźniki błędu maksymalnego");
xlabel('N');
ylabel('K');
zlabel('err_max');

%% Zadanie 4
clear;
N = 40;
T = 99;
vectorsigm=logspace(-4,-1,20);
err = zeros(31,1);
result = zeros(T,31);
result_srednia = zeros(31,numel(vectorsigm));
K_opt = zeros(numel(vectorsigm),1);

for p=1:numel(vectorsigm)
    sigma=vectorsigm(p);
    for j =1:T
        S = sigma*randn(N,1);
        for i=5:35
             result(j,i-4) = error_zaburzone(N,i,S);
        end
    end
    [~,index] = min(mean(result));
    K_opt(p,1) = index + 4;
end

figure(5);
plot(log10(vectorsigm), K_opt, "*b");
title("Zadanie 4");
ylabel('K_opt');
xlabel('log10(Sigma)');
grid on;

%% Zadanie 4 - przyklady
clear;
N = 40;
x = linspace(-1,1,1000);
vectorsigm=logspace(-4,-1,20)
sigma = vectorsigm(20);
S = sigma*randn(N,1);
[Q_1000, Q_N] = aprox_zaburzenie(40,30,S);

figure(6);
zad1(N);
hold on;
plot(x,Q_1000, "b");
plot(linspace(-1,1,N),Q_N, "ok")
xlabel('x');
ylabel('f(x)');
title(sprintf('Zadanie 4 - przykład'));
grid on;
legend('Orygnalna funkcja f(x)', 'Aproksymacja', "Location", "northwest")
hold off;
%% Funkcje
%Rysowanie wykresu podanej funkcji
function [f]= zad1(N)
    x = linspace(-1,1,1000);
    f = (2/3).*cos(pi.*(x-1/3)).*exp(x);
    f_x = linspace(-1,1,N);
    f_k = (2/3).*cos(pi.*(f_x-1/3)).*exp(f_x);

    plot(x ,f, 'r');
    hold on;
    %plot(f_x, f_k, 'bo');
    hold off;
end

%Funkcja obliczajaca wartosci wspolczynnikow p
function [M,Y,P] = rysowanie(N,K)
    f_x = linspace(-1,1,N);
    x_k = linspace(-1,1,K);
    f_k = (2/3).*cos(pi.*(f_x-1/3)).*exp(f_x);

    Y = zeros(N,1);
    M = zeros(N,K);
    for i=1:K
        Gx = (4.*exp(-16.*(f_x-x_k(i)).*(f_x-x_k(i))))./sqrt(2.*pi);
        for j=1:N
            M(j,i)=Gx(j);
            Y(j,1)=f_k(j);
        end
    end
    P = (M'*M)\(M'*Y);
end

%Funkcja aproksymujaca wykres, tj. mnozenie funkcji przez wspolczynniki p
function [Q_1000,Q_N] = aprox(N,K)
    [M_N,~,P] = rysowanie(N,K);
    [M_1000,~,~] = rysowanie(1000,K);
    Q_1000 = zeros(1000,1);
    F_1000 = zeros(1000,K);
    Q_N = zeros(N,1);
    F_N = zeros(N,K);
    for i=1:1000
        for j=1:K
            F_1000(i,j) = M_1000(i,j) * P(j);
            Q_1000(i) = Q_1000(i) + F_1000(i,j);
        end
    end
    for i=1:N
        for j=1:K
            F_N(i,j) = M_N(i,j) * P(j);
            Q_N(i) = Q_N(i) + F_N(i,j);
        end
    end
end

%Funkcja liczaca wskaźnik błędu średniokwadratowego
function err = error(N, K)
    x = linspace(-1,1,1000);
    f = (2/3).*cos(pi.*(x-1/3)).*exp(x);
    [Q_1000, ~] = aprox(N,K);

    R = zeros(1000,1);
    sum = 0;
    arg = 0;
    for i=1:1000
        R(i,1) = abs(Q_1000(i,1) - f(i));
        sum = sum + abs(R(i,1)).^2;
        arg = arg + abs(f(i)).^2;
    end
    err = sqrt(sum)./arg;
end

%Funkcja liczaca wskaźnik błędu maksymalnego
function err_max = error_max(N, K)
    x = linspace(-1,1,1000);
    f = (2/3).*cos(pi.*(x-1/3)).*exp(x);
    [Q_1000, ~] = aprox(N,K);

    R = zeros(1000,1);
    for i=1:1000
        R(i,1) = (Q_1000(i,1) - f(i));
    end
    err_max = max(R)./max(f);
end

%Funkcja obliczajaca wartosci wspolczynnikow p
function [M,Y,P] = zaburzenie(N,K,S)
    f_x = linspace(-1,1,N);
    x_k = linspace(-1,1,K);
    f_k = (2/3).*cos(pi.*(f_x-1/3)).*exp(f_x);


    Y = zeros(N,1);
    M = zeros(N,K);
    for i=1:K
        Gx = (4.*exp(-16.*(f_x-x_k(i)).*(f_x-x_k(i))))./sqrt(2.*pi);
        for j=1:N
            M(j,i)=Gx(j);
            Y(j,1)=f_k(j) + S(j,1);
        end
    end
    P = (M'*M)\(M'*Y);
end

%Funkcja aproksymujaca wykres, tj. mnozenie funkcji przez wspolczynniki p
function [Q_1000,Q_N] = aprox_zaburzenie(N,K,S)
    [M_N,~,P] = zaburzenie(N,K,S);
    [M_1000,~,~] = rysowanie(1000,K);
    Q_1000 = zeros(1000,1);
    F_1000 = zeros(1000,K);
    Q_N = zeros(N,1);
    F_N = zeros(N,K);
    for i=1:1000
        for j=1:K
            F_1000(i,j) = M_1000(i,j) * P(j);
            Q_1000(i) = Q_1000(i) + F_1000(i,j);
        end
    end
    for i=1:N
        for j=1:K
            F_N(i,j) = M_N(i,j) * P(j);
            Q_N(i) = Q_N(i) + F_N(i,j);
        end
    end
end

%Funkcja liczaca wskaźnik błędu średniokwadratowego
function err = error_zaburzone(N, K, S)
    x = linspace(-1,1,1000);
    f = (2/3).*cos(pi.*(x-1/3)).*exp(x);
    [Q_1000, ~] = aprox_zaburzenie(N,K, S);

    R = zeros(1000,1);
    sum = 0;
    arg = 0;
    for i=1:1000
        R(i,1) = abs(Q_1000(i,1) - f(i));
        sum = sum + abs(R(i,1)).^2;
        arg = arg + abs(f(i)).^2;
    end
    err = sqrt(sum)./arg;
end