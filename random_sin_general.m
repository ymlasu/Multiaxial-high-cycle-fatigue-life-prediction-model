function [F_life, sig_spec, tau_spec, max_angle, critical_angle1] = random_sin_general(sig, tau, sig_mean, tau_mean, phy, s, Faxial, Baxial, sig_ult)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  This code is developed by Haoyang Wei and Yongming Liu.       %%%%%
%%%%%%  All content can be found in the published article.            %%%%%
%%%%%%  If you use this code or any of the included functions for     %%%%%
%%%%%%  scientific purpose please respect the effort and cite the     %%%%%
%%%%%%  paper which name is shown below.                              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Title: Multiaxial high-cycle fatigue life prediction under random spectrum loadings %%%

%%% Web: https://doi.org/10.1016/j.ijfatigue.2019.105462 %%%

%%% Author: Haoyang Wei, Yongming Liu, and other co-author %%%
%%% Arizona State University, AZ %%%

%% Material properties -- The detail can be found at Table.1
cosine_2alfa = (-2 + sqrt(4-4*(1/s^2-3)*(5-1/s^2-4*s^2))) / 2 / (5-1/s^2-4*s^2);
alfa = acos(cosine_2alfa)/2;
beta = sqrt(cosine_2alfa^2*s^2 + (1-cosine_2alfa^2));
eta = 1.5;

%% Tensile and torsion loading spectrum  
% Build constant loading history
fs = 60;
angle_interval = 6;
angles = [0: 2*pi/360*angle_interval: 2*pi];
sig_peak = sig*sin(angles + 0.5*pi) + sig_mean;
tau_peak = tau*sin(angles+phy/180*pi + 0.5*pi) + tau_mean;

lX = length(angles)-1;
sig_spec = -diff(sig_peak)/2.*cos(pi*(0:1/fs:1-1/fs)') + (sig_peak(1:lX)+sig_peak(2:lX+1))/2;
sig_spec = [sig_spec(:);sig_peak(end)];
tau_spec = -diff(tau_peak)/2.*cos(pi*(0:1/fs:1-1/fs)') + (tau_peak(1:lX)+tau_peak(2:lX+1))/2;
tau_spec = [tau_spec(:);tau_peak(end)];

%% Find critical plane angle (maximum tensile damage plane +- alfa) 
% Liu-Mahadevan critical plane concept can be found at chapter 2.1
angle_step = 1;
damage_max = 0;
damage_max_angle = 0;
for angle = 0: 2*pi/360*angle_step: 2*pi
    sig_angle = sig_spec/2 + sig_spec/2*cos(2*angle) + tau_spec*sin(2*angle);
    c_sig_angle = rainflow(sig_angle,fs);
    damage_sig_angle = c_sig_angle(:,1)' * (1./(10.^(10.^((log10(c_sig_angle(:,2)/2)-log10(Faxial))./(Baxial)))));
    if damage_sig_angle > damage_max_angle
        damage_max_angle = damage_sig_angle;
        max_angle = angle;
    end
end
critical_angle1 = max_angle + alfa;
critical_angle2 = max_angle - alfa;
sig_max_ang = sig_spec/2 + sig_spec/2*cos(2*max_angle) + tau_spec*sin(2*max_angle);

%%  tensile and torsion damage on critical plane, equivalent damage
% The detail can be found at chapter 2.2 and equation start from Eq.1
sig_critical1 = sig_spec/2 + sig_spec/2*cos(2*critical_angle1) + tau_spec*sin(2*critical_angle1);
tau_critical1 = -sig_spec/2*sin(2*critical_angle1) + tau_spec*cos(2*critical_angle1);

c_sig_critical1 = rainflow(sig_critical1);
c_tau_critical1 = rainflow(tau_critical1);

damage_sig_critical1 = c_sig_critical1(:,1)' * c_sig_critical1(:,2)/2;
damage_tau_critical1 = c_tau_critical1(:,1)' * c_tau_critical1(:,2)/2/s;

damage_sig_critical1 = (1./(10.^(10.^((log10(damage_sig_critical1 * (1 + eta * sum(sig_max_ang)/size(sig_max_ang,1) / sig_ult))-log10(Faxial))./(Baxial)))));
damage_tau_critical1 = (1./(10.^(10.^((log10(damage_tau_critical1 * (1 + eta * sum(sig_max_ang)/size(sig_max_ang,1) / sig_ult))-log10(Faxial))./(Baxial)))));

sig_eqv1 = 10 ^ (Baxial * log10(0- log10(damage_sig_critical1)) + log10(Faxial));
tau_eqv1 = 10 ^ (Baxial * log10(0 - log10(damage_tau_critical1)) + log10(Faxial));

stress_eqv1 = sqrt(sig_eqv1^2 + tau_eqv1^2)/beta;
damage_eqv1 = 10^0 * (1./(10.^(10.^((log10(stress_eqv1)-log10(Faxial))./(Baxial)))));

sig_critical2 = sig_spec/2 + sig_spec/2*cos(2*critical_angle2) + tau_spec*sin(2*critical_angle2);
tau_critical2 = -sig_spec/2*sin(2*critical_angle2) + tau_spec*cos(2*critical_angle2);

c_sig_critical2 = rainflow(sig_critical2);
c_tau_critical2 = rainflow(tau_critical2);

damage_sig_critical2 = c_sig_critical2(:,1)' * c_sig_critical2(:,2)/2;
damage_tau_critical2 = c_tau_critical2(:,1)' * c_tau_critical2(:,2)/2/s;

damage_sig_critical2 = (1./(10.^(10.^((log10(damage_sig_critical2 * (1 + eta * sum(sig_max_ang)/size(sig_max_ang,1) / sig_ult))-log10(Faxial))./(Baxial)))));
damage_tau_critical2 = (1./(10.^(10.^((log10(damage_tau_critical2 * (1 + eta * sum(sig_max_ang)/size(sig_max_ang,1) / sig_ult))-log10(Faxial))./(Baxial)))));

sig_eqv2 = 10 ^ (Baxial * log10(0 - log10(damage_sig_critical2)) + log10(Faxial));
tau_eqv2 = 10 ^ (Baxial * log10(0 - log10(damage_tau_critical2)) + log10(Faxial));

stress_eqv2 = sqrt(sig_eqv2^2 + tau_eqv2^2)/beta;
damage_eqv2 = 10^0 * (1./(10.^(10.^((log10(stress_eqv2)-log10(Faxial))./(Baxial)))));

if damage_eqv1 > damage_eqv2
    damage_eqv = damage_eqv1;
else
    damage_eqv = damage_eqv2;
end

%% Fatigue life based on equivalent damage
F_life = log10(1/damage_eqv);
end
