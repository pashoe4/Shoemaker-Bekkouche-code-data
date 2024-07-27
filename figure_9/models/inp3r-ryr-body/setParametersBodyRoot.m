function out = setParametersBodyRoot(infile,outfile,varargin)
%%%% setParametersBodyRoot.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a cell body model,
%   with attached processes & process roots
% A setParametersProc script MUST BE RUN PRIOR TO RUNNING THIS SCRIPT:
%   cell body parameters utilize basic cellular process parameters
%   NOTE: Laplacian matrix fo processes MUST BE 2nd ORDER!
% User adjusts user-settable parameters in script
% When run, script prompts for unique output filename for .mat file
% Cell thickness is VARIABLE in this formulaton

%% Parse input parameters
p = inputParser;
addOptional(p,'volRb',0.15);% volume ratio, cytosol/ER 
addOptional(p,'kEb',4);% ER fold factor
%addOptional(p,'fCacb',0.2);%Ca diffusion reduction parameter. Fraction to reduce cyto. DCa due to intracellular crowding
%addOptional(p,'fCaeb',0.4);%Ca diffusion reduction parameter. Fraction to reduce ER DCa due to intracellular crowding
addOptional(p,'RrIb',10);% InP3R density, units um^-2
%addOptional(p,'fInP3b',0.33);%InP3 diffusion reduction parameter. Fraction to reduce DInP3 due to intracellular crowding
addOptional(p,'fCoVar',1);%Factor used to co-vary fInP3b, fCaer and fCacr
addOptional(p,'np',6);%number of processes
parse(p,varargin{:});
volRb=p.Results.volRb;
kEb=p.Results.kEb;
%fCacb=p.Results.fCacb;
%fCaeb=p.Results.fCaeb;
RrIb=p.Results.RrIb;
%fInP3b=p.Results.fInP3b;
fCoVar=p.Results.fCoVar;
np=p.Results.np;

%infile = input('Enter name of parameter file for cellular process: ','s');

%outfile = input('Specify name of output file: ','s');

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
tbo = 3; % outer edge (units um)
tbi = 6; % adjacent to nucleus (units um)

% number of cellular processes
%np = 6;

% root length (note: process dx must be even divisor)
lr = 2; % units um
% root half-width (increase over process radius rd at junction)
wr = 0.6; % units um

% ER fold factor, cell body
%kEb = 5;
% ER fold factor, process roots
kEr = kEb/2;

% Ca TRANSPORT/REACTION PARAMETERS

% Ca diffusion reduction parameters
fCacr = 0.2*fCoVar; % fraction to reduce radial cyto. DCa, intracell. crowding
fCaer = 0.4*fCoVar; % fraction to reduce radial ER DCa due to intracell. crowding

% InP3R density
%RrIb = 10; % units um^-2

% InP3 diffusion reduction parameter
fInP3r = 0.33*fCoVar; % fraction to reduce radial DInP3 due to intracell. crowding

%% Computation-related parameters

% MORPHOLOGICAL / GRID PARAMETERS

% # vol elements in circumferential direction
nc = 48; % SHOULD BE INTEGER MULTIPLE OF np! ****
ac = 2*pi()/nc; % angle subtended by vol elements (units radians)

nct = nc/np; % # circumferential elements from junction to junction ****
ncp=ceil(2*(rd+wr)/(ac*rb)); % # circmf elts subtended by root @ jnctn ****

nr = 6;  % # vol elements in radial direction ****
nrc = 1:nr-1; % indices of central radial vol elements ****
drb = (rb-rn)/nr; % radial grid spacing

tb = (tbi:-(tbi-tbo)/(nr-1):tbo); % total elt. thickness as fnctn of radius
nt = 10; % # vol elements in z (thickness) direction
dtb = tb./nt; % vol element thicknesses

ngr = lr/dx; % # process root compartments
dwr = wr/ngr; % increment in width per compartment

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
    rxb = [rxb, rxb(end)+drb];
end
rc = rxb(1:nr) + 0.5*drb; % centerline radii

% PM area of sidewalls
apsw = 2*pi()*rb*tbo / (nt*nc);

% PM areas of element faces
for k=1:nr
    apmf(k) = pi()*(rxb(k+1)^2 - rxb(k)^2) / nc;
end

% element ER areas
aet = kEb*2*pi()*mean(rc)*mean(tb); % total ER membrane area
aem = rc.*tb./sum(rc.*tb) * aet / (nt*nc); % ER area per element

% element volumes
for k=1:nr
    volb(k) = ( rxb(k+1)^2-rxb(k)^2 ) * tb(k) / (nt*nc); % total volume
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
RBW = R0./drb^2;
RBJ = R1./drb^2;
RBJb = R1b./drb^2;

% constants for (2nd order) z (thickness)-direction diffusion computations
T0 = zeros(nt,nt); 
% for sealed end at bottom
T0(1,2) = 1;
T0(1,1) = -1;
% for middle elements
for k=2:nt-1
    T0(k,k-1) = 1;
    T0(k,k+1) = 1;
    T0(k,k) =  -2;
end
% for sealed end at top
T0(end,end-1) = 1;
T0(end,end) = -1;
TB = cell(nr,1); % cell array to hold scaled diffusion coeff. for z-dir.
for k=1:nr
    TB{k} = T0./dtb(k)^2; % scaling by (z-direction length)^-2
end

% get boundary condition coefficients
RB = robinBC(kD, nr, drb, 0, 1); % kth-order Neumann BC constants
% row matrix to compute right side grad BC, with BV coeff. separated out
JB = full(RB(end,:));
JBb = JB(end);
JB = JB(2:end-1);

% process root parameters

% cross-sectional areas of radial elements
for k=1:3
    axcr(k) = pi()*( rxp(k+1)*(rxp(k+1)+wr) - rxp(k)*(rxp(k)+wr) );
end

% area fractional weights for process root cytoplasmic elements
fxp = axcr./sum(axcr);

% expand process to include root grid points

kpr = 1:ng; % indices for process proper
krt = ng+1:ng+ngr; % indices for root

% average diffusion scaling for root
wrk = rd + (0:ngr-1)*dwr;
brl = (rd+wrk -   4*drp) ./ (rd+wrk+0.5*dwr-4*drp);
brr = (rd+wrk+dwr-4*drp) ./ (rd+wrk+0.5*dwr-4*drp);

%   L1 = matrix for last process grid point and ngr root grid points
L1 = lap(2,ngr+2,dx); % sparse Laplacian matrix, order 2, increment dx
for k=4:ngr+3 % rescale coefficients relating to root
    L1(k,k-1) = brl(k-3)*L1(k,k-1);
    L1(k,k+1) = brr(k-3)*L1(k,k+1);
    L1(k,k) = -( L1(k,k-1)+L1(k,k+1) );
end
L1 = L1(3:end,:); % get rid of first two rows
L1 = [ zeros(ngr+2,ng-2),L1 ]; % pad with zeros on left

% expand LP to include L1
LP = LP(1:ng,:);  % get rid of last two rows
LP = [ LP,zeros(ng,ngr) ]; % pad with zeros on right
LP = [ LP; L1 ]; % assemble new LP matrix

% expand boundary value coefficient vectors
BR = [ zeros(1,ngr), BR ];
BL = [ BL, zeros(1,ngr) ];
JP = [ zeros(1,ngr), JP ];

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%

% Calcium molar membrane leakage rates
% rate to cytosol from extracellular space
%   for non-sidewall volume elements
klxbf = rlx * atvpmf * kvol; % units uM.s^-1 ****
%   for volume elements with sidewalls
klxbw = rlx * atvpmw * kvol; % units uM.s^-1 ****
% rate to cytosol from body ER
kleb = rle * atver * kvol; % units s^-1 ****
% rate to cytosol from root ER (rescaled from process parameter)
kler = kEr * kle; % units s^-1 ****

% Ca diffusion coefficients
% in cytosol
DCacr = DCa*fCacr; % units um^2.s^-1 ****
% in ER
DCaer = DCa*fCaer; % units um^2.s^-1 ****

%%% Calcium pump parameters %%%

% Calcium NCX pump constants (values differ for ea volume element)
%   for face areas (note these can differ with radius)
NpNf = round(RpN*apmf);
kpNbf = NpNf*IpN ./volbc * kvol; % units uM.s^-1 ****
%   for sidewalls
NpNw = round(RpN*apsw*nt); % mult.then div.by nt: 'de-quantize' indiv elts
kpNbw = (NpNw*IpN/nt) ./volbc(end) * kvol;  % units uM.s^-1 ****

% Calcium PMCA pump constant
%   for face areas
kpPbf = RpP*IpP * atvpmf * kvol; % units uM.s^-1 ****
%   for sidewalls
kpPbw = RpP*IpP * atvpmw * kvol; % units uM.s^-1 ****

% Calcium SERCA pump constants
kpSb =  RpS*IpS * atver * kvol;  % units uM.s^-1 ****
% for roots (rescaled from process parameter)
kpSr = kEr * kpS; % units s^-1 ****

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3r = DInP3*fInP3r;  % units um^2.s^-1 ****

% Ca molar flow rate through InP3R channels
%   vector: values differ for elts of ea radius
NrIb = round(RrIb*aem*nt); % mult.then div.by nt: 'de-quantize' indiv elts
krIb = (NrIb*IrI/nt) ./volbc * kvol; % units s^-1. ****
% for roots (rescaled from process parameter)
krIr = kEr * krI;

% rate constant, [Ca]-dependent production of InP3
%   for non-sidewall volume elements
kmIfbf = rmIf * atvpmf * kvol; % multiplies [Ca]; units s^-1 ****
%   for volume elements with sidewalls
kmIfbw = rmIf * atvpmw * kvol; % multiplies [Ca]; units s^-1 ****

%%% Ryanodine receptor parameters %%%
% excluded at present

%% Write output

clear k dr rxb apsw apmf aet aem volb atvpmf atvpmw atver...
      L0 L1 R0 R1 R1b RB axcr wrk brl brr kvol

save(outfile);
