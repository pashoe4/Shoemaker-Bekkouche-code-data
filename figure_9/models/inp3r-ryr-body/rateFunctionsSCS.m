%%%% rateFunctions.m %%%%
% Script to compute & store receptor rate constants and functions
% Includes simplified Siekmann-Cao-Sneyd InP3 model, B&Q RyR model

% Notation, output parameters:
%   ph used for Ca-dependent rate functions
%   r  used for rate constants
%   f = forward; b = backward (reaction direction relative to reference)

%clear

%% Siekman-Cao-Sneyd inositol trisphosphate receptor model:
% Notes:
% (pos. and neg. digits used in paper for forward & backward, resp.);
% REFERENCE DIRECTION: forward <=> *away* from R state
%   1st row: paper notation, SKS (digits are subscripts);
%   2nd row: my notation below
%   q12   q21  q23  q32   q24   q42   q26  q62  q45  q54
%   rI12 rI21 rI23 rI32  phI24 phI42 rI26 rI62 rI45 rI54
% forward (f) => park to drive; back (b) => drive to park
%
% Original model structure
%  
%                  OI6
%                  ||
%              rI26||rI62
%           rI12   ||  rI23
%       CI1 ====== CI2 ====== CI3
%           rI21   ||  rI32
%                  ||        'Drive mode'
% ----------- phI42||phI24 --------------
%                  ||        'Park mode'
%           rI45   ||
%       OI5 ====== CI4 
%           rI54
%
% Modified model structure:
%
%           rI26       rI21
%       OI6 ====== CI2 ===== CI1
%           rI62   ||  rI12
%                  ||        'Drive mode'
% ----------- phI42||phI24 --------------
%                  ||        'Park mode'
%           rI45   ||  phI43
%       OI5 ====== CI4 ===== CI3
%           rI54       rI34
%

% receptor kinetic rate constants
rI12 =  1240;   % units s^-1
rI21 =    88;   % units s^-1
rI26 = 10500;   % units s^-1
rI62 =  4010;   % units s^-1
%rI45 =    11;   % units s^-1 % original value, SCS
rI45 =   2.2;   % units s^-1 % reduced value
rI54 =  3330;   % units s^-1
rI34 =  0.05;   % units s^-1

% parameter associated with mode transition rate phI43
rI43 =  1000;  % units s^-1.uM^-2

% parameters associated with mode transition rates phI42 & phI24
    % (rate constants, corresponding to V's of Cao et al.
rI42  =   50;  % units s^-1
rI24  =  400;  % units s^-1


% gating parameters for Ca dependence of phRs
km    =   0.08;  % units uM
km6pw =  km.^6;  % units uM^6
kh    =    3.5;  % units uM
kh6pw =  kh.^6;  % units uM^6

    % gating parameters for InP3 dependence of phRs
kbp   =    0.5;  % units uM
% kbp4pw = kbp.^4; % units uM^4
kbp3pw = kbp.^3; % units uM^3

% NOTE: rate constant related to gating functions is set in
%       setParameters.m, to allow experimental variation


%% Breit & Queisser ryanidine receptor model
% Notes: 'pos' and 'neg'
% REFERENCE DIRECTION: forward <=> *away* from C1 state
%   1st row: paper notation, B&Q (+/- are superscripts, digits subscripts);
%   2nd row: my intermediate notation below
%   k+a k-a k+b k-b k+c k-c kab kbb kcb kcf  O1 O2
%   kaf kab kbf kbb kcf kcb rab rbb rcb rcf  O3 O4
%
%             OR4
%             ||
%        phR34||rR43
%      phR13  ||     rR32
% CR1 ======  OR3  ======  CR2
%      rR31          rR23
%    

% B-Q fundamental channel constants
rR31 = 28.8;     % units s^-1
rR43 = 385.9;    % units s^-1
rR23 = 0.1;      % units s^-1
rR32 = 1.75;     % units s^-1

% parameters associated with calcium-dependent rates phR13 & phR34
kR13 = 1500;     % units uM^-4s^-1
kR34 = 1500;     % units uM^-3s^-1

% parameters for saturation of phR13 & phR34
rsat34 = 10*(rR31+rR43+rR23+rR32); % saturation value for phR13
Ccsat = (rsat34./kR34).^(1/3); % [Ca] value at which this occurs;
rsat13 = kR13.*Ccsat.^4; % corresponding saturation value for phR34

%% Save computed data

save InP3R_rates.mat rI12 rI21 rI26 rI62 rI45 rI54 rI34 rI43 rI42 rI24...
                     km6pw kh6pw kbp3pw
                 
save RyR_rates.mat rR31 rR43 rR23 rR32 kR13 kR34 rsat13 rsat34

%save('InP3R_rates','rR31', 'rR43', 'rR23', 'rR32', 'kR13', 'kR34', 'rsat13', 'rsat34');