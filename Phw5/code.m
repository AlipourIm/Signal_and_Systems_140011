%%
% Section 1: Laplace Transform

% Problem 1:
%%
syms time

% part a:
f_1_t = time*heaviside(time-1);   % signal for part a
laplace(f_1_t)

%%

% part b:
f_2_t = sin(time)*exp(-4*time) * heaviside(time);     % signal for part b
laplace(f_2_t)

%%

% part c:
f_3_t = 2*time*cos (3*time)*heaviside(time);     % signal for part c
laplace(f_3_t)

%%
% Problem 2:

syms s

% part a:
F_1_s = (1/(s*(s+1)))*exp(-3*s);   % signal for part a
ilaplace(F_1_s)

%%

% part b:
F_2_s = 4/(s*(s^2 +4));    % signal for part b
ilaplace(F_2_s)

%%

% part c:
F_3_s = 1/(s^2 + 3 * s + 1);   % signal for part c
ilaplace(F_3_s)

%% 
% Problem 3:

G_s =(8)/(s^2 +s + 4);
impulse_response = ilaplace(G_s)
step_response = ilaplace(G_s*(1/s))

%% 
% Part a: impulse response
time = -1:0.001:10;
y_1 = subs(impulse_response,time) .* (time>=0);
plot_the_figure(time,y_1,"impulse response h(time)")

%% 
% Part a: step response
time2 = -1:0.001:10;
y_2 = subs(step_response,time2) .* (time2>=0);
plot_the_figure(time2,y_2,"step response H(s)")

%% 
% Part b:   Bode plot of H(s)
H_s = tf((8),[1,1,4]);
bode(H_s);
%% 
% Problem 4:

% For a = 4 answer is:

G_s =(2*s+1)/(s^2 + 4*s + 7);
step_response = ilaplace(G_s*(1/s));
time = -1:0.001:10;
y_2 = subs(step_response,time) .* (time>=0);
plot_the_figure(time,y_2,"step response with a = 4");

clear time
syms time

step_response =  sym(1/7) - (exp((-2*time))*(cos(sqrt(sym(3))*time) - 4*sqrt(sym(3))*sin(sqrt(sym(3))*time)))/7;
disp("Limit at infinity = ");
limit(step_response,time,Inf);
time= 0:0.001:2;
y = subs(step_response,time) .* (time>=0);
[value,index_of] = max(y);

% display the results;
disp("Maximum is at = ");
disp(time(index_of))
disp("Maximum value is :")
disp(double(value))
%% 

% Problem 4:

% For a = 6

G_s =(2*s+1)/(s^2 + 6*s + 7);
step_response = ilaplace(G_s*(1/s));
time = -1:0.001:10;
y_2 = subs(step_response,time) .* (time>=0);
plot_the_figure(time,y_2,"step response with a = 6");

clear time
syms time

step_response =  sym(step_response);
disp("Limit at infinity = ");
limit(step_response,time,Inf)
time= 0:0.001:2;
y = subs(step_response,time) .* (time>=0);
[value,index_of] = max(y);

% Display the results:
disp("Maximum is at = ");
disp(time(index_of));
disp("Maximum value is :");
disp(double(value));

%% Section 2: Z transform
% 
% Problem 1:




%% Problem 1:
% Part a:

syms n
x1 = (heaviside(n+2)-heaviside(n-1))*(-n+1)-(heaviside(n-1) - heaviside(n-3))*(n+1);
n1 = -5:5;
step_figure_maker(n1,subs(x1,n1),"plot");
x1 = (heaviside(n)-heaviside(n-3))*(-n+3)+(heaviside(n-3) - heaviside(n-5))*(-n+1);
n1 = -5:10;
step_figure_maker(n1,subs(x1,n1),"plot");
syms z
z_transform_plan_plotter([3,2,1,-2,-3],[1,0,0]);

%% Problem 1: part b

% Part b:
x2 = 0.8^n * heaviside(n-2);
t2 = 0:15;
step_figure_maker(t2,subs(x2,t2),"plot");
z_transform_plan_plotter([5/4,0],[5/4,-1,0,0]);

%% Problem 1: part c

% Part c:
x3 = 2^n* cos(0.4*pi*n) * heaviside(n);
t3 = -2:6;
step_figure_maker(t3,subs(x3,t3),"plot");
H = (z-2*cos(0.4*pi))/(z^2-4*cos(0.4*pi)*z + 4)
z_transform_plan_plotter([1,-2*cos(0.4*pi)],[1,-4*cos(0.4*pi),4]);
[r,p,k] = residuez([1,-2*cos(0.4*pi)],[1,-4*cos(0.4*pi),4]);
%% Problem 2:
% 
% 
% Part a

syms z
num1 = [1,-1];
dom1 = [1,-1,0.5];
H1 = (z-1)/(z^2 - z + 0.5)
z_transform_plan_plotter(num1,dom1);

num2 = [1,0];
dom2 = [2, - sqrt(3),0.5];
H2 = (z)/(z^2 - sqrt(3)*z + 0.5)
z_transform_plan_plotter(num2,dom2);
%% 
% Part b:

[r1,p1,k1] = residuez(num1,dom1);
[r2,p2,k2] = residuez(num2,dom2);
%% 
% Part c:

old_parameter = sympref('HeavisideAtOrigin',1);
h1 = iztrans(H1)
h2 = iztrans(H2)
%% 
% Helper Functions: 

function plot_the_figure(x,y,titlf)
figure;
set(gcf,'position',[0,0,1800,900]);
plot(x,y);
title(titlf);
shg;
end

function z_transform_plan_plotter(num,dom)
H = tf(num,dom);
figure;
set(gcf,'position',[0,0,1800,900]);
pzplot(H);
h = findobj(gca, 'type', 'line');
set(h, 'markersize', 20)
text(real(roots(num)) - 0.1, imag(roots(num)) + 0.1, 'Zero')
text(real(roots(dom)) - 0.1, imag(roots(dom)) + 0.1, 'Pole')
axis equal
figure;
set(gcf,'position',[0,0,1800,900]);
zplane(num,dom);
grid on;
end

function step_figure_maker(x,y,titlf)
figure;
set(gcf,'position',[0,0,1800,900]);
stem(x,y);
title(titlf);
end