% Building a reference field model from Swarm data

clear all
a = 6371.2; % earth radius [km]
rad = pi/180; % from deg to rad
%% (i) Load and extract the data
load('data.mat','Blambda','Br','Btheta','F','lambda','r','t','theta');
%% (ii) and (iii) Build a design matrix and find a least-squares solution
[Gr, Gtheta, Glambda] = design_SHA(r/a, theta*rad, lambda*rad, 13);
m = [Gr; Gtheta; Glambda]\[Br; Btheta; Blambda];  % least-squares solution
g01 = m(1); 
rate = (abs(g01-(-29446)))/4; % rate of change
fprintf('The design matrix has 10114 rows ans 195 columns\n')
fprintf('It changed %1.4f nT on average per year during\n',rate)
%% (iv) Plot histograms of the residuals between your model predictions and the Swarm data
Br_pred = Gr*m;
Btheta_pred = Gtheta*m;
Blambda_pred = Glambda*m;

figure
subplot(3,2,1)
plot(theta, Br, '.r')
title('B_r from data'), xlabel('θ [deg]'), ylabel('B_r [nT]')
subplot(3,2,2)
plot(theta, Br_pred, '.b')
title('B_r predicted'), xlabel('θ [deg]'), ylabel('B_r [nT]')
subplot(3,2,3)
plot(theta, Btheta, '.r')
title('B_θ from data'), xlabel('θ [deg]'), ylabel('B_r [nT]')
subplot(3,2,4)
plot(theta, Btheta_pred, '.b')
title('B_θ predicted'), xlabel('θ [deg]'), ylabel('B_r [nT]')
subplot(3,2,5)
plot(theta, Blambda, '.r')
title('B_λ from data'), xlabel('θ [deg]'), ylabel('B_r [nT]')
subplot(3,2,6)
plot(theta, Blambda_pred, '.b')
title('B_λ predicted'), xlabel('θ [deg]'), ylabel('B_r [nT]')

% residuals

DeltaBr = Br-Br_pred;
DeltaBtheta = Btheta-Btheta_pred;
DeltaBlambda = Blambda-Blambda_pred;
F_pred = sqrt((Br_pred.^2)+(Btheta_pred.^2)+(Blambda_pred.^2));
DeltaF = F-F_pred;

figure
subplot(2,2,1)
histogram(DeltaBr)
title('ΔB_r'), xlabel('Residuals'), ylabel('Relative frequency')
subplot(2,2,2)
histogram(DeltaBtheta)
title('ΔB_θ'), xlabel('Residuals'), ylabel('Relative frequency')
subplot(2,2,3)
histogram(DeltaBlambda)
title('ΔB_λ'), xlabel('Residuals'), ylabel('Relative frequency')
subplot(2,2,4)
histogram(DeltaF)
title('ΔF'), xlabel('Residuals'), ylabel('Relative frequency')

figure
subplot(2,2,1)
plot(theta, DeltaBr, 'ob')
xlabel('θ [deg]'), ylabel('ΔB_r [nT]')
subplot(2,2,2)
plot(theta, DeltaBtheta, 'ob')
xlabel('θ [deg]'), ylabel('ΔB_θ [nT]')
subplot(2,2,3)
plot(theta, DeltaBlambda, 'ob')
xlabel('θ [deg]'), ylabel('ΔB_λ [nT]')
subplot(2,2,4)
plot(theta, DeltaF, 'ob')
xlabel('θ [deg]'), ylabel('ΔF [nT]')

% rsm
rsm_r = sqrt(1/length(Br)*(sum(DeltaBr.^2)));
rsm_theta = sqrt(1/length(Btheta)*(sum(DeltaBtheta.^2)));
rsm_lambda = sqrt(1/length(Blambda)*(sum(DeltaBlambda.^2)));
fprintf('The RSM misfit for each of the three components: %2.4f, %2.4f and %2.4f\n%',rsm_r,rsm_theta,rsm_lambda)
%% (v) Plot maps of the intensity, magnetic declination and magnetic inclination
[theta, lambda] = meshgrid(1:179, -180:2:180);
theta = reshape(theta, [], 1);
lambda = reshape(lambda, [], 1);
[Gr, Gtheta, Glambda] = design_SHA(1, theta*rad, lambda*rad, 13);

Br_new = Gr*m;
Btheta_new = Gtheta*m;
Blambda_new = Glambda*m;
X = -Btheta_new;
Y = Blambda_new;
Z = -Br_new;
F_new = sqrt((X.^2)+(Y.^2)+(Z.^2));
H = sqrt((X.^2)+(Y.^2));
D = rad2deg(atan(Y./X));
I = rad2deg(atan(Z./H));

load('world');
figure
subplot(3,1,1)
scatter(lambda, theta, [], F_new, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet;
colorbar;
title('F [nT]');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')

subplot(3,1,2)
scatter(lambda, theta, [], D, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet; 
colorbar;
title('Declination');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')

subplot(3,1,3)
scatter(lambda, theta, [], I, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet; 
colorbar;
title('Inclination');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')
%% (vi) Magnetic declination D in Copenhagen
[Gr, Gtheta, Glambda] = design_SHA(1, 34.32*rad, 12.56*rad, 13);
Br_CPH = Gr*m;
Btheta_CPH = Gtheta*m;
Blambda_CPH = Glambda*m;
X_CPH = -Btheta_CPH;
Y_CPH = Blambda_CPH;
D_CPH = rad2deg(atan(Y_CPH./X_CPH));
fprintf('The magnetic declination D in Copenhagen in sep 2018 is %1.4f degrees\n',D_CPH)
%% (vii) Model with truncation degree N = 50 vs IGRF-12
load('data.mat','Blambda','Br','Btheta','F','lambda','r','t','theta');
[Gr, Gtheta, Glambda] = design_SHA(r/a, theta*rad, lambda*rad, 50);
m_50 = [Gr; Gtheta; Glambda]\[Br; Btheta; Blambda];
[theta, lambda] = meshgrid(1:179, -180:2:180);
theta = reshape(theta, [], 1);
lambda = reshape(lambda, [], 1);
[Gr, Gtheta, Glambda] = compute_sha_matrices(1, theta*rad, lambda*rad, 50);

Br_50 = Gr*m_50;
Btheta_50 = Gtheta*m_50;
Blambda_50 = Glambda*m_50;
X_50 = -Btheta_50;
Y_50 = Blambda_50;
Z_50 = -Br_50;
F_50 = sqrt((X_50.^2)+(Y_50.^2)+(Z_50.^2));

[theta, lambda] = meshgrid(1:179, -180:2:180);
theta = reshape(theta, [], 1);
lambda = reshape(lambda, [], 1);
[Gr, Gtheta, Glambda] = design_SHA(1, theta*rad, lambda*rad, 13);
m_igrf12 = load('igrf12.txt');
Br_igrf12 = Gr*m_igrf12;
Btheta_igrf12 = Gtheta*m_igrf12;
Blambda_igrf12 = Glambda*m_igrf12;

X_igrf12 = -Btheta_igrf12;
Y_igrf12 = Blambda_igrf12;
Z_igrf12 = -Br_igrf12;
F_igrf12 = sqrt((X_igrf12.^2)+(Y_igrf12.^2)+(Z_igrf12.^2));

f = F_igrf12-F_50;

figure
subplot(3,1,1)
scatter(lambda, theta, [], F_50, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet;
colorbar;
title('N=50 [nT]');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')

subplot(3,1,2)
scatter(lambda, theta, [], F_igrf12, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet;
colorbar;
title('igrf12 [nT]');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')

subplot(3,1,3)
scatter(lambda, theta, [], f, 'filled');
axis([-180 180 0 180])
set(gca, 'YDir', 'reverse');
colormap jet;
colorbar;
title('differences [nT]');
xlabel('λ [deg]');
ylabel('θ [deg]');
hold on
plot(coastline.X,90-coastline.Y,'k')
%% (viii) Co-estimate a dipolar external field along with your internal field model
load('data.mat','Blambda','Br','Btheta','F','lambda','r','t','theta');
[Gr, Gtheta, Glambda] = design_SHA(r/a, theta*rad, lambda*rad, 13);
[Gr_ext, Gtheta_ext, Glambda_ext] = design_SHA(r/a, theta*rad, lambda*rad, 1, 'ext');
% m_ext = [Gr_ext; Gtheta_ext; Glambda_ext]\[Br; Btheta; Blambda];
Gr_tot = [Gr, Gr_ext];
Gtheta_tot = [Gtheta, Gtheta_ext];
Glambda_tot = [Glambda, Glambda_ext];
G_tot = [Gr_tot; Gtheta_tot; Glambda_tot];
m_tot = [G_tot]\[Br; Btheta; Blambda];