function [F_life, damage_eqv, critical_angle1, max_angle] = spectrum_prediction(sig_spec,tau_spec,s,Faxial,Baxial, sig_ult)

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

%% Find critical plane angle (maximum tensile damage plane +- alfa)
% Liu-Mahadevan critical plane concept can be found at chapter 2.1
angle_step = 1;
damage_max = 0;
damage_max_angle = 0;
for angle = 0: 2*pi/360*angle_step: 2*pi
    sig_angle = sig_spec/2 + sig_spec/2*cos(2*angle) + tau_spec*sin(2*angle);
    c_sig_angle = rainflow(sig_angle);
    damage_sig_angle = c_sig_angle(:,1)' * (1./(10.^(10.^((log10(c_sig_angle(:,2)/2)-log10(Faxial))./(Baxial)))));
    
    if damage_sig_angle > damage_max_angle
        damage_max_angle = damage_sig_angle;
        max_angle = angle;
    end
end
critical_angle1 = max_angle + alfa;
critical_angle2 = max_angle - alfa;
sig_max_angle = sig_spec/2 + sig_spec/2*cos(2*max_angle) + tau_spec*sin(2*max_angle);

%% Tensile and torsion damage on critical plane, equivalent damage
% The detail can be found at chapter 2.2 and equation start from Eq.1
sig_critical1 = sig_spec/2 + sig_spec/2*cos(2*critical_angle1) + tau_spec*sin(2*critical_angle1);
tau_critical1 = -sig_spec/2*sin(2*critical_angle1) + tau_spec*cos(2*critical_angle1);

% Rainflow acounting algorithm
c_sig_critical1 = rainflow(sig_critical1);
c_tau_critical1 = rainflow(tau_critical1);


for i = 1:size(c_sig_critical1)
    mean_ten1(i,1) = (sum(sig_max_angle(c_sig_critical1(i,4) : c_sig_critical1(i,5))) - 0.5 * (sig_max_angle(c_sig_critical1(i,5)) + sig_max_angle(c_sig_critical1(i,4)))) / (c_sig_critical1(i,5) - c_sig_critical1(i,4));
end

for i = 1:size(c_tau_critical1,1)
    mean_shear1(i,1) = (sum(sig_max_angle(c_tau_critical1(i,4) : c_tau_critical1(i,5))) - 0.5 * (sig_max_angle(c_tau_critical1(i,5)) + sig_max_angle(c_tau_critical1(i,4)))) / (c_tau_critical1(i,5) - c_tau_critical1(i,4));
end

damage_sig_critical1 = c_sig_critical1(:,1)' * (1./(10.^(10.^((log10(c_sig_critical1(:,2)/2 .* (1 + eta * mean_ten1 / sig_ult))-log10(Faxial))./(Baxial)))));
damage_tau_critical1 = c_tau_critical1(:,1)' * (1./(10.^(10.^((log10(c_tau_critical1(:,2)/2/s .* (1 + eta * mean_shear1 / sig_ult))-log10(Faxial))./(Baxial)))));

sig_eqv1 = Faxial * 10^(Baxial * log10(log10(size(c_sig_critical1,1) / damage_sig_critical1)));
tau_eqv1 = Faxial * 10^(Baxial * log10(log10(size(c_sig_critical1,1) / damage_tau_critical1)));

stress_eqv1 = sqrt(sig_eqv1^2 + tau_eqv1^2)/beta;
damage_eqv1 = size(c_sig_critical1,1) * (1./(10.^(10.^((log10(stress_eqv1)-log10(Faxial))./(Baxial)))));

sig_critical2 = sig_spec/2 + sig_spec/2*cos(2*critical_angle2) + tau_spec*sin(2*critical_angle2);
tau_critical2 = -sig_spec/2*sin(2*critical_angle2) + tau_spec*cos(2*critical_angle2);

c_sig_critical2 = rainflow(sig_critical2);
c_tau_critical2 = rainflow(tau_critical2);

for i = 1:size(c_sig_critical2)
    mean_ten2(i,1) = (sum(sig_max_angle(c_sig_critical2(i,4) : c_sig_critical2(i,5))) - 0.5 * (sig_max_angle(c_sig_critical2(i,5)) + sig_max_angle(c_sig_critical2(i,4)))) / (c_sig_critical2(i,5) - c_sig_critical2(i,4));
end

for i = 1:size(c_tau_critical2,1)
    mean_shear2(i,1) = (sum(sig_max_angle(c_tau_critical2(i,4) : c_tau_critical2(i,5))) - 0.5 * (sig_max_angle(c_tau_critical2(i,5)) + sig_max_angle(c_tau_critical2(i,4)))) / (c_tau_critical2(i,5) - c_tau_critical2(i,4));
end

damage_sig_critical2 = c_sig_critical2(:,1)' * (1./(10.^(10.^((log10(c_sig_critical2(:,2)/2 .* (1 + eta * mean_ten2 / sig_ult))-log10(Faxial))./(Baxial)))));
damage_tau_critical2 = c_tau_critical2(:,1)' * (1./(10.^(10.^((log10(c_tau_critical2(:,2)/2/s .* (1 + eta * mean_shear2 / sig_ult))-log10(Faxial))./(Baxial)))));

sig_eqv2 = Faxial * 10^(Baxial * log10(log10(size(c_sig_critical2,1) / damage_sig_critical2)));
tau_eqv2 = Faxial * 10^(Baxial * log10(log10(size(c_sig_critical2,1) / damage_tau_critical2)));

stress_eqv2 = sqrt(sig_eqv2^2 + tau_eqv2^2)/beta;
damage_eqv2 = size(c_sig_critical2,1) * (1./(10.^(10.^((log10(stress_eqv2)-log10(Faxial))./(Baxial)))));

if damage_eqv1 > damage_eqv2
    damage_eqv = damage_eqv1;
else
    damage_eqv = damage_eqv2;
end

%% Fatigue life based on equivalent damage
F_life = log10(1/damage_eqv*length(sig_spec));
end
