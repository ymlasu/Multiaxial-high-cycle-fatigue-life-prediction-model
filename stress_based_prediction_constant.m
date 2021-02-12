%% Multiaxial high-cycle fatigue life prediction for constant loading

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
load('constant_loading.mat'); % File from the experimental data.
Load = constant_loading;
F_life = [];
angle = [];

% Tesion and torsion fitting parameters
a_ten = 2176.2; b_ten = -1.351; 
a_tor = 1495.3; b_tor = -1.357; 
s = (a_tor * 4^b_tor)/(a_ten * 4^b_ten);

% Yield strength
sig_y = 503;

% Fatigue prediction 
for i = 1:length(Load)
    [F_life(i,1),~,~,~,angle(i,1)] = random_sin_general(Load(i,1),Load(i,2),Load(i,3), Load(i,4), Load(i,5),s, a_ten, b_ten, sig_y);
end

% Predicted fatigue life
F_life