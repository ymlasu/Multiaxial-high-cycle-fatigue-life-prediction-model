%% Multiaxial high-cycle fatigue life prediction for random loading

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

clc, clear all
% Data input and pre-setting
load('Proportional_Felix.mat'); % File from the experimental data.

Load = Load;
F_life = [];
angle = [];

% Tesion and torsion fitting parameters
a_ten = 2176.2; b_ten = -1.351;
a_tor = 1495.3; b_tor = -1.357;
s = (a_tor * 4^b_tor)/(a_ten * 4^b_ten);

% Yield strength
sig_y = 503;

% Fatigue prediction 
[F_life,damage,criticle_angle,failure_angle] = spectrum_prediction(Load(:,1), Load(:,2), s, a_ten, b_ten, sig_y);

% Predicted fatigue life
F_life

