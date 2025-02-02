%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\bo.bekkouche\Documents\code\NMDAwave\calciumWaveModels-9.5\extra_analysis\ParamDepend_Graphs.xlsx
%    Worksheet: Ampl noRyR
%
% Auto-generated by MATLAB on 13-May-2023 12:28:56

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 24);

% Specify sheet and range
opts.Sheet = "Ampl noRyR";
opts.DataRange = "A2:X10";

% Specify column names and types
opts.VariableNames = ["rlxCeCtl", "Ampl", "rle", "Ampl1", "fCac", "Ampl2", "fCae", "Ampl3", "fInP3", "Ampl4", "Rp_PM", "Ampl5", "RpS", "Ampl6", "kbCf", "Ampl7", "rmIf", "Ampl8", "kiPI", "Ampl9", "rIg", "Ampl10", "Ce0", "Ampl11"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
ampl_noryr = readtable("C:\Users\bo.bekkouche\Documents\code\NMDAwave\calciumWaveModels-9.5\extra_analysis\ParamDepend_Graphs.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts