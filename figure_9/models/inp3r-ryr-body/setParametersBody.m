function out = setParametersBody(infile,outfile,varargin)
%%%% setParametersBody.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a cell body model,
%   with attached processes
% A setParametersProc script MUST BE RUN PRIOR TO RUNNING THIS SCRIPT:
%   cell body parameters utilize basic cellular process parameters
% User adjusts user-settable parameters in script
% When run, script prompts for unique output filename for .mat file
% Cell thickness is VARIABLE in this formulaton

%clear;

%% Parse input parameters
p = inputParser;
addOptional(p,'volRb',0.15);% volume ratio, cytosol/ER 
addOptional(p,'kEb',5);% ER fold factor
addOptional(p,'fCacb',0.2);%Ca diffusion reduction parameter. Fraction to reduce cyto. DCa due to intracellular crowding
addOptional(p,'fCaeb',0.4);%Ca diffusion reduction parameter. Fraction to reduce ER DCa due to intracellular crowding
addOptional(p,'RrIb',10);% InP3R density, units um^-2
addOptional(p,'fInP3b',0.33);%InP3 diffusion reduction parameter. Fraction to reduce DInP3 due to intracellular crowding
parse(p,varargin{:});
volRb=p.Results.volRb;
kEb=p.Results.kEb;
fCacb=p.Results.fCacb;
fCaeb=p.Results.fCaeb;
RrIb=p.Results.RrIb;
fInP3b=p.Results.fInP3b;
%infile = 'parameters_proc';%input('Enter name of parameter file for cellular process: ','s');

%outfile = 'parameters_body';%input('Specify name of output file: ','s');

%% Read in parameters for mathing cellular process model

load(infile);

%% User-settable (at present) model parameter values

% MORPHOLOGICAL PARAMETERS

% volume ratio, cytosol/ER 
%volRb = 0.15;

% cell body radius:
%   allows process w/ r=1.5um to intersect 4x7.5deg circumf. body units
rb = 6; % units um
% nuclear radius
rn = 2.4; % units um

% cell body thickness
tbo = 3; % units um
tbi = 6; % units um

% number of cellular processes
%np = 8;
np = 6;

% ER fold factor
%kEb = 5;

% Ca TRANSPORT/REACTION PARAMETERS

% Ca diffusion reduction parameters
%fCacb = 0.2; % fraction to reduce cyto. DCa due to intracellular crowding
%fCaeb = 0.4; % fraction to reduce ER DCa due to intracellular crowding

% InP3R density
%RrIb = 10; % units um^-2

% InP3 diffusion reduction parameter
%fInP3b = 0.33; % fraction to reduce DInP3 due to intracellular crowding

%% Computation-related parameters

% GRID PARAMETERS

% # vol elements in circumferential direction
nc = 48; % SHOULD BE INTEGER MULTIPLE OF np! ****
ac = 2*pi()/nc; % angle subtended by vol elements (units radians)

nct = nc/np; % # circumferential elements from junction to junction ****
ncp=ceil(2*rd/(ac*rb)); % # circumf elts subtended by process @ jnctn ****

nr = 6;  % # vol elements in radial direction ****
nrc = 1:nr-1; % indices of central radial vol elements ****
dr = (rb-rn)/nr; % radial grid spacing
tb = (tbi:-(tbi-tbo)/(nr-1):tbo); % element thicknesses

% circumferential indices of junction & non-junction vol elements 
ksw = [];
kjn = [];
k0 = 0;
for k=1:np % loop over processes
    kjn = [ kjn, k0 + (1:ncp) ];
    ksw = [ ksw, k0 + (ncp+1:nct) ];
    k0 = k0 + nct;
end

% boundary radii for cytoplasmic elements
rxb = rn;
for k=1:nr
    rxb = [rxb, rxb(end)+dr];
end
rc = rxb(1:nr) + 0.5*dr; % centerline radii

% PM area of sidewalls
apsw = 2*pi()*rb*tbo / nc;

% PM areas of element faces
for k=1:nr
    apmf(k) = 2*pi()*(rxb(k+1)^2 - rxb(k)^2) / nc;
end

% element ER areas
aet = kEb*2*pi()*mean(rc)*mean(tb); % total ER membrane area
aem = rc.*tb./sum(rc.*tb) * aet / nc; % ER area per element

% element volumes
for k=1:nr
    volb(k) = ( rxb(k+1)^2-rxb(k)^2 ) * tb(k) / nc; % total volume
end
volbe = volRb/(volRb+1) * volb; % ER volumes
volbc = 1/(volRb+1) * volb; % cyto volumes

% element PM area/volume ratios
atvpmf = apmf./volbc; % for faces
atvpmw = apsw/volbc(end); % for sidewalls

% element ER membrane area / cyto volume ratio (same at all radii)
atver = aem(1)/volbc(1);

% Constants for circumferential diffusion computations (DERIVED WITH MOLE) 
kD = 2; % order of Laplacian computation (presumed preset)
L1=lap(kD,2*kD+1,1); % minimal-sized sparse Laplacian matrix, order kD
L1 = full(L1);
L0 = L1(kD+2,:); % middle row of L1
L0 = L0(3:end-2); % vector of (2*kD-1) non-zero Laplacian coefficients

% scaling by (circumferential 'length')^-2
LB = cell(nr,1); % cell array to hold scaled Laplacian coefficients
for k=1:nr
    dc = sin(ac)*rc(k);
    LB{k} = L0./dc^2;
end

% constants for (2nd order) radial diffusion computations
R0 = zeros(nr,nr); 
% for sealed end at inner radius
R0(1,2) = 2 * tb(2)/tb(1) * rxb(2)/(rxb(1)+rxb(2));
R0(1,1) = -R0(1,2);
% for middle elements
for k=2:nr-1
    R0(k,k-1) = 2 * rxb(k)/(rxb(k)+rxb(k+1));
    R0(k,k+1) = 2 * tb(k+1)/tb(k) * rxb(k+1)/(rxb(k)+rxb(k+1));
    R0(k,k) = -( R0(k,k-1) + R0(k,k+1) );
end
% constants for process junction elements (same except for last row)
R1 = R0;
% last row for sealed end at outer radius
R0(end,end-1) = 2 * rxb(end-1)/(rxb(end-1)+rxb(end));
R0(end,end) = -R0(end,end-1);
% last row for process junction at outer radius
R1(end,end-1) = 2.6667*rxb(nr)/(rxb(nr)+rxb(nr+1));
R1b = 5.3333*rxb(nr+1)/(rxb(nr)+rxb(nr+1));
R1(end,end) = -( R1(end,end-1) + R1b );
% scaling by (radial length)^-2
RBW = R0./dr^2;
RBJ = R1./dr^2;
RBJb = R1b./dr^2;

% get boundary condition coefficients
RB = robinBC(kD, nr, dr, 0, 1); % kth-order Neumann BC constants
% row matrix to compute right side grad BC, with BV coeff. separated out
JB = full(RB(end,:));
JBb = JB(end);
JB = JB(2:end-1);

% area fractional weights for process root cytoplasmic elements
fxp = axcp./sum(axcp);

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%

% Calcium molar membrane leakage rates
% rate to cytosol from extracellular space
%   for non-sidewall volume elements
klxbf = rlx * atvpmf * kvol; % units uM.s^-1 ****
%   for volume elements with sidewalls
klxbw = rlx * atvpmw * kvol; % units uM.s^-1 ****
% rate to cytosol from ER
kleb = rle * atver * kvol; % units s^-1 ****

% Ca diffusion coefficients
% in cytosol
DCacb = DCa*fCacb; % units um^2.s^-1 ****
% in ER
DCaeb = DCa*fCaeb; % units um^2.s^-1 ****

%%% Calcium pump parameters %%%

% Calcium NCX pump constants (values differ for ea volume element)
%   for face areas (note these can differ with radius)
NpNf = round(RpN*apmf);
kpNbf = NpNf*IpN ./volbc * kvol; % units uM.s^-1 ****
%   for sidewalls
NpNw = round(RpN*apsw);
kpNbw = NpNw*IpN ./volbc(end) * kvol;  % units uM.s^-1 ****

% Calcium PMCA pump constant
%   for face areas
kpPbf = RpP*IpP * atvpmf * kvol; % units uM.s^-1 ****
%   for sidewalls
kpPbw = RpP*IpP * atvpmw * kvol; % units uM.s^-1 ****

% Calcium SERCA pump constants
kpSb =  RpS*IpS * atver * kvol;  % units uM.s^-1 ****

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3b = DInP3*fInP3b;  % units um^2.s^-1 ****

% Ca molar flow rate through InP3R channels
%   vector: values differ for elts of ea radius
NrIb = round(RrIb*aem);
krIb = (NrIb*IrI) ./volbc * kvol; % units s^-1. ****

% rate constant, [Ca]-dependent production of InP3
%   for non-sidewall volume elements
kmIfbf = rmIf * atvpmf * kvol; % multiplies [Ca]; units s^-1 ****
%   for volume elements with sidewalls
kmIfbw = rmIf * atvpmw * kvol; % multiplies [Ca]; units s^-1 ****

%%% Ryanodine receptor parameters %%%
% excluded at present

%% Write output

clear k dr rxb apsw apmf aet aem volb atvpmf atvpmw atver...
      L0 L1 R0 R1 R1b RB kvol p varargin

save(outfile);
