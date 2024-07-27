% Simulates and saves data for Figure 9. Also generates plots for each
% variable.

clear all
close all

addpath(genpath('utils'))
addpath(genpath('models'))
addpath(genpath('experiments'))
addpath(genpath('analyses'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text_label= '_test6'; % Name label for the results folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate body root and vary parameters (cell body edge to edge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative change to variable. This is the x-axis of Figure 9.
fVar=0.25:0.25:2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate, analyse and plot each variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='RrIb';results_name=['exp_BodyRootCaDynamics_RrIb_tuned' text_label];interval=10*fVar;var_location='setParametersBodyRoot';var_str_leg='R_{rI}(Ce ctl)'
% Path to tuning values:
kerval_path=[pwd '\models\inp3r-ryr-body\ker_val_map.mat'];%Won't use if ''
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kerval_path='';
var_str='volRb';results_name=['exp_BodyRootCaDynamics_volRb' text_label];interval=0.15*fVar;var_location='setParametersBodyRoot';var_str_leg=var_str;
%exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='kEb';results_name=['exp_BodyRootCaDynamics_kEb' text_label];interval=4*fVar;var_location='setParametersBodyRoot';var_str_leg=var_str;
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='RrIb';results_name=['exp_BodyRootCaDynamics_RrIb' text_label];interval=10*fVar;var_location='setParametersBodyRoot';var_str_leg=var_str;
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='kbCf';results_name=['exp_BodyRootCaDynamics_kbCf' text_label];interval=1.9*fVar;var_location='setParametersProc';var_str_leg=var_str;
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='rmIf';results_name=['exp_BodyRootCaDynamics_rmIf' text_label];interval=1.4E-15*fVar;var_location='setParametersProc';var_str_leg=var_str;
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='rIg';results_name=['exp_BodyRootCaDynamics_rIg' text_label];interval=70*fVar;var_location='setParametersProc';var_str_leg=var_str;
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='fCoVar';results_name=['exp_BodyRootCaDynamics_fCoVar' text_label];interval=1*fVar;var_location='setParametersBodyRoot';var_str_leg='f_{CoVar}';
exp_BodyRootCaDynamics(results_name,var_str,interval,var_location,kerval_path);% Simulate
ana_BodyRootCaDynamics_body(results_name,var_str,interval,kerval_path,var_str_leg);% Analyse and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the variations above together (subfigures of Figure 9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str='';var_str_list={'volRb','kEb','RrIb','RrIb','kbCf','rmIf','rIg','fCoVar'};
var_str_leg={'vol_{Rb}','k_{Eb}','R_{rI}','R_{rI}(Ce ctl)','k_{bCf}','r_{mIf}','r_{Ig}','f_{CoVar}'};
intervals={0.15*fVar, 4*fVar, 10*fVar, 10*fVar, 1.9*fVar,1.4E-15*fVar,70*fVar,fVar};results_name=text_label;%
tuned=[4,9999];%index of the tuned variable (RrIb) in list above
default_i=4;%index of fVar==1

% Element index on morphology where measurements are made
elem_nr = 9;
ana_BodyRootCaDynamics_multi_body(results_name,var_str_list,intervals,fVar,tuned,elem_nr,default_i,var_str_leg);

elem_nr = 17;
ana_BodyRootCaDynamics_multi_body(results_name,var_str_list,intervals,fVar,tuned,elem_nr,default_i,var_str_leg);

elem_nr = 25;%Default at bottom
ana_BodyRootCaDynamics_multi_body(results_name,var_str_list,intervals,fVar,tuned,elem_nr,default_i,var_str_leg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The plots have been generated next to this file. In addition, for each
% simulation, data and snapshot plots have been generated in the results
% subfolders.