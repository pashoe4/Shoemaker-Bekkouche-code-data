%%%% ProcCaDynamicsTaper3D_ranRyR.m %%%%
% Compute equilibrium state of a compartment with ER + RYRs, followed by
% time-domain response of a tapered cylindrical process with ER + RYRs,
% with 3-D diffusion.
% 'Excitation' takes form of pulsed Ca & InP3 fluxes at distal end.
% InP3 model is Siekmann-Cao-Sneyd model with re-defined CI3 state;
% RyR model is from Breit & Quiesser / Keizer & Levine.
% RyRs may be discretely distributed along dendrite, < 1/compartment,
% randomly distributed around the circumference of the ER
%
% User is prompted for a file containing desired parameter values.

clear;

%dt = 1E-5;   % time step for temporal integration (units s)
dt = 5E-7;   % small time step for temporal integration (units s)

iplot = 500; % plot interval in units of dt
tout = iplot*dt; % plot interval

%% Load parameters

% load model parameters
mpfile = input('Specify name of model parameter file: ','s');
load(mpfile); % tapered process model parameters

% load rate constants / function tables for InP3 receptors
load('InP3R_rates.mat');

% load rate constants / function tables for Ryanodine receptors
load('RyR_rates.mat');

%% Set distal Ca & InP3 influx parameters for activated processes

% influx rate
Cain = 200; % units uM.s^-1
InP3in = 200; % units uM.s^-1

BVCain = BCc * Cain;
BVInP3in = BInP3 * InP3in;

% influx pulse duration
tin = 0.02; % units s

%% Scale rate constants by dt for purposes of temporal integration

%%% Universal calcium parameters %%%

% Calcium molar leakage rates
DClx = dt*klx; % (note extracellular leakage assumed constant)
kle = dt*kle;

% Cytosol Ca diffusion coefficient
DCac = dt*DCac;
% ER Ca diffusion coefficient
DCae = dt*DCae;

% calcium buffering rates
kbCf = dt*kbCf;
kbCb = dt*kbCb;

%%% Calcium pump parameters %%%

% pump clearance rate constants
kpN = dt*kpN;
kpP = dt*kpP;
kpS = dt*kpS;

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = dt*DInP3;

% Calcium molar influx rate for open InP3R channels
krI = dt*krI;

% InP3 generation rate
kmIf = dt*kmIf;
% InP3 degradation rate
kmIb = dt*kmIb;

% InP3 receptor state rates
rI12 = dt*rI12;  
rI21 = dt*rI21;
rI26 = dt*rI26;
rI62 = dt*rI62;
rI45 = dt*rI45;
rI54 = dt*rI54;
rI34 = dt*rI34;

rI43 = dt*rI43;
rI42 = dt*rI42;
rI24 = dt*rI24;

% rate parameter for InP3 receptor gating
rIg = rIg*dt;

% PKC activation rate
kiPf = dt*kiPf;
% PKC deactivation rate
kiPb = dt*kiPb;

%%% Ryanodine receptor parameters %%%

% Calcium molar influx rate for open RyR channels
krR = dt*krR;

% RyR receptor state rates
rR23 = dt*rR23;
rR32 = dt*rR32;
rR31 = dt*rR31;
rR43 = dt*rR43;

kR13 = dt*kR13;
kR34 = dt*kR34;

rsat13 = dt*rsat13;
rsat34 = dt*rsat34;

% Input flux increments
Cain = dt*Cain;
InP3in = dt*InP3in;

%% Take means of parameters for equilibrium calculation

DClxX = mean(DClx);
kleX = mean(kle);
kpNX = mean(kpN);
kpPX = mean(kpP);
kpSX = mean(kpS);
kmIfX = mean(kmIf);
krIX = mean(krI);
krRX = mean(krR);
MrRX = nRyR/ngt; % fraction of compartments with RyRs

%% Set up initial conditions / initial states for equilibrium calculation

% Small calcium concentrations used to ensure convergence to equilibrium
CcX = 0.01; % cytosolic Ca; units uM
CeX = 200; % ER Ca; units uM

CaCalBX = CcX*kbCf*CalB0/kbCb; % bound calcium value

% approximate quiescent InP3-related states
InP3X = kmIfX*CcX/kmIb; % cytosolic InP3;
CI1X = 0;
CI2X = 0;
CI3X = 0;
CI4X = 1;
OI5X = 0;
OI6X = 0;
PKCvX = 0;
% gating variables
mItX = 0;
hItX = 1;

% allocate InP3R state increments
DI21 = 0;
DI26 = 0;
DI42 = 0;
DI45 = 0;
DI43 = 0;

% approximate quiescent RyR states
CR1X = 1;
CR2X = 0;
OR3X = 0;
OR4X = 0;

% allocate RyR state increments
DR13 = 0;
DR32 = 0;
DR34 = 0;

% allocate InP3R calcium increment
DCI = 0;

% allocate RyR calcium increments
DCR = 0;

%% Relaxation to observe (unexcited) equilibrium states, or instability

% output variables to examine history
Cct=[]; Cet = []; time = [];

it = 0;

DClx0 = volRoc*DClxX; % rescaled external leakage for single compartment
toutr = 20*tout; % downsample outputs

tend = 100; % end time
for t = 0:dt:tend % loop on time
    
    % Compute state increments
    
    % Ca leakage from ER
    DCle = volRic*kleX.*(CeX-CcX);
    
    % Ca buffering in cytosol
    DCfb = kbCf.*CcX.*(CalB0-CaCalBX);
    DCbf = kbCb.*CaCalBX;
    
    % Calcium pumps
    
        % NCX pumps
    DCN = volRoc*kpNX.*CcX./(kpNC+ CcX);
    
    % PMCA pumps
    Cc2pw =  CcX.* CcX;
    DCP = volRoc*kpPX.*Cc2pw./(kpPC2pw+Cc2pw);
    
    % SERCA pumps
    DCS = volRic*kpSX.* CcX./((kpSC+ CcX).*CeX);
    
    % InP3R and related states
    
    % Gating function dynamics
    Cc6pw = CcX.^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3X.^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mIt = mI;
%     hIt = hI;

    % single-pole delay
        % (done here since these play role in InPrR state updates)
    mItX = mItX + rIg.*( mI - mItX );
    hItX = hItX + rIg.*( hI - hItX );

    % InP3 states increments
    phI43 = rI43.*(1 - hItX);
    phGate = mIp.*mItX.*hItX;
    phI42 = rI42.*(    phGate);
    phI24 = rI24.*(1 - phGate);
    DI21 = rI21.*CI2X - rI12.*CI1X;
    DI26 = rI26.*CI2X - rI62.*OI6X;
    DI42 = phI42.*CI4X - phI24.*CI2X;
    DI45 = rI45.*CI4X - rI54.*OI5X;
    DI43 = phI43.*CI4X - rI34.*CI3X;
    
    % probability of open state for InP3Rs
    OpenI = OI5X + OI6X;
    % Ca increment through open InP3R channels
    DCI = volRic*krIX*OpenI.*( CeX - CcX );
    
    % cytosolic InP3 increment
    DInP3c = (volRoc*kmIfX* CcX)./(1+kiPI*PKCvX) - kmIb*InP3X;
    
    % PKCv update
    DPKCv = kiPf* CcX.*(PKC0-PKCvX) - kiPb*PKCvX;
    
    % RyR states
    
    rR13 = kR13.*CcX.^4;
    rR34 = kR34.*CcX.^3;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);

    % RyR states increments
    DR13 = phR13.*CR1X - rR31.*OR3X;
    DR32 = rR32.*OR3X - rR23.*CR2X;
    DR34 = phR34.*OR3X - rR43.*OR4X;
    
    % probability of open state for RyR
    OpenR = OR3X + OR4X;
    
    % Ca increment through open RyR channels
        % note scaling by 1/nc for circumferential divisions
    DCR = volRic*krRX*OpenR.*( CeX - CcX )./(nc*MrRX);
    
    % Apply increments to states
    
    % sum of ER calcium increments
    DCe = ( DCS - DCR-DCI - DCle )/volR;
    % update Ce
    CeX = CeX + DCe;
    
    % sum of cytosolic calcium increment
    DCc = DCbf-DCfb - DCN-DCP-DCS + DCI+DCR + DClx0+DCle;
    % update Cc
    CcX = CcX + DCc;
    
    % InP3 concentration increment
    InP3X = InP3X+DInP3c;
    
    % Calbindin state
    CaCalBX = CaCalBX + DCfb - DCbf;
    
    % InP3R and related states
    CI1X = CI1X + DI21;
    CI2X = CI2X + DI42 - DI21 - DI26;
    CI3X = CI3X + DI43;
    CI4X = CI4X - DI45 - DI42 - DI43;
    OI5X = OI5X + DI45;
    OI6X = OI6X + DI26;
    
    PKCvX = PKCvX + DPKCv;
    
    % RyR states
    CR1X = CR1X - DR13;
    CR2X = CR2X + DR32;
    OR3X = OR3X + DR13 - DR32 - DR34;
    OR4X = OR4X + DR34;
    
    if ~mod(t,toutr)
        time = [time;t];
        Cct = [Cct,CcX];
        Cet = [Cet,CeX];
    end
    
end

figure(11)
plot(time, Cct,'b-')
grid on
xlabel('Time (s)');
ylabel('Concentration (uM)');
set(gca,'fontsize',15);
figure(12)
plot(time, Cet,'-','color', [0,0.3,0.6])
grid on
xlabel('Time (s)');
ylabel('Concentration (uM)');
set(gca,'fontsize',15);

Cc_Cmprt = CcX;
Ce_Cmprt = CeX;
InP3_Cmprt = InP3X;

%% Set up initial conditions / initial states for dendrite

% load arrays for dendrite with approximate quiescent parameter values

Cc = CcX*ones(ngt+2,3,nc); % cytosolic Ca
Ce = CeX*ones(ngt+2,1); % ER Ca
InP3 = InP3X*ones(ngt+2,3,nc); % cytosolic InP3;
CaCalB = CaCalBX*ones(ngt,3,nc); % bound calcium value

% quiescent InP3R states
CI1 = CI1X*ones(ngt,1,nc);
CI2 = CI2X*ones(ngt,1,nc);
CI3 = CI3X*ones(ngt,1,nc);
CI4 = CI4X*ones(ngt,1,nc);
OI5 = OI5X*ones(ngt,1,nc);
OI6 = OI6X*ones(ngt,1,nc);
% gating variables
mIt = mItX*ones(ngt,1,nc);
hIt = hItX*ones(ngt,1,nc);

% activated PKC
PKCv = PKCvX*ones(ngt,1,nc);

% quiescent RyR  states
CR1 = CR1X*ones(nRyR,1);
CR2 = CR2X*ones(nRyR,1);
OR3 = OR3X*ones(nRyR,1);
OR4 = OR4X*ones(nRyR,1);

% pre-allocate arrays for 3D increments
DCcRx = zeros(3,nc,ngt);
DInP3Rx = zeros(3,nc,ngt);
DCcCx = zeros(nc,3,ngt);
DInP3Cx = zeros(nc,3,ngt);
DCcA = zeros(ngt+2,3,nc);
DInP3A = zeros(ngt+2,3,nc);

% allocate array for ER leak increments
DCle = zeros(ngt,1,nc);

% allocate arrays for pump calcium increments
DCS = zeros(ngt,1,nc);
DCP = zeros(ngt,1,nc);
DCN = zeros(ngt,1,nc);

% allocate arrays for InP3R state increments
DI21 = zeros(ngt,1,nc);
DI26 = zeros(ngt,1,nc);
DI42 = zeros(ngt,1,nc);
DI45 = zeros(ngt,1,nc);
DI43 = zeros(ngt,1,nc);

% allocate arrays for RyR state increments
DR13 = zeros(nRyR,1);
DR32 = zeros(nRyR,1);
DR34 = zeros(nRyR,1);

% allocate arrays for receptor calcium increments
DCI = zeros(ngt,1,nc);
DCR = zeros(ngt,1,nc);

%% Set up for time loop

% grids for plots
xgrid1 = [0.5*dx:dx:(ngt-0.5)*dx]';
xlims = [-1,ldt+1];
xtic = [0:ldt/10:ldt];

Cct = {};  % cell array for calcium history
it = 0; % counter for indexing history
ngstart = ng0;  % index into CcM at start of taper
ngend = ng0+ng+1; % index into CcM at end of taper

%% Loop on time

tend = 0.75; % end time

for t = 0:dt:tend % loop on time
    
    Cc(Cc<0) = 0; % added to prevent numerical undershoot in [Ca]
    InP3(InP3<0) = 0; % added to prevent numerical undershoot in [InP3]

%% internal Cc, Ce, & InP3R for computing states @ grid points; aux. arrays

    % internal Cc, InP3 values
    CcM = Cc(2:end-1,:,:);
    InP3M = InP3(2:end-1,:,:);

    % permute 2nd & 3rd (radial & circumf.) indices of CcM, InP3M
    CcMx = permute(CcM,[1 3 2]); % temp array (to be reused later)
    InP3Mx = permute(InP3M,[1 3 2]); % temp array (to be reused later)
    % 2D arrays for inner radial elements
    CcMy = CcMx(:,:,1); 
    InP3My = InP3Mx(:,:,1);
    % averages of permuted arrays around circumference
    CcMave = pagetranspose(mean(pagetranspose(CcMx)));
    InP3Mave = pagetranspose(mean(pagetranspose(InP3Mx)));
    % permute indices back to standard order
    CcMave = permute(CcMave,[1 3 2]);
    InP3Mave = permute(InP3Mave,[1 3 2]);
    % save as 2D arrays
    CcMave = CcMave(:,:,1);
    InP3Mave = InP3Mave(:,:,1);
    
        % Internal Ce values
    CeM = Ce(2:end-1,1);
    
    
%% graphing and archiving
    
     if ~mod(t,tout)

% archive CcMave data
        it = it+1;
        Cct{it} = CcMave(ngstart:ngend,:);
        
% graphing
        % set up figure for cytosolic Ca
        figure(1)
        plot(xgrid1,CcMave(:,1),'g.-');
        hold on
        grid on
        plot(xgrid1,CcMave(:,2),'.-','color',[0,0.75,0.75]);
        plot(xgrid1,CcMave(:,3),'b.-');
        plot(xgrid1,InP3Mave(:,1),'r-')
        xlabel('Distance (um)');
        ylabel('Concentration (uM)');
        set(gca,'fontsize',15);
        xlim(xlims);
        xticks(xtic);
        ylim([0 20]);
        hold off
        
        % set up figure for ER Ca
        figure(2)
        plot(xgrid1,CeM(:,1),'.-','color', [0.3,0.3,0.8]);
        hold on
        grid on
        xlabel('Distance (uM)');
        ylabel('Concentration (uM)');
        set(gca,'fontsize',15);
        xlim(xlims);
        xticks(xtic);
        ylim([0 500]);
        hold off

     end

%% Compute increments for all states %%

    %% Ca leakage from ER

    DCle = kle.*(CeM-CcM(:,1,:));

    %% Ca buffering in cytosol
    
    DCfb = kbCf.*CcM.*(CalB0-CaCalB);
    DCbf = kbCb.*CaCalB;
 
    %% Calcium pumps
    
    % NCX pumps
    DCN = kpN.*CcM(:,3,:)./(kpNC+CcM(:,3,:));
    
    % PMCA pumps
    Cc2pw = CcM(:,3,:).*CcM(:,3,:);
    DCP = kpP.*Cc2pw./(kpPC2pw+Cc2pw);

    % SERCA pumps
    DCS = kpS.*CcM(:,1,:)./((kpSC+CcM(:,1,:)).*CeM);

    %% InP3R and related states
        
    % Gating function dynamics
    Cc6pw = CcM(:,1,:).^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3M(:,1,:).^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mIt = mI;
%     hIt = hI;

    % single-pole delay
        % (done here since these play role in InPrR state updates)
    mIt = mIt + rIg.*( mI - mIt );
    hIt = hIt + rIg.*( hI - hIt );
    
    % InP3R states increments
    phI43 = rI43.*(1 - hIt);
    phGate = mIp.*mIt.*hIt;
    phI42 = rI42.*(    phGate);
    phI24 = rI24.*(1 - phGate);
    DI21 = rI21.*CI2 - rI12.*CI1;
    DI26 = rI26.*CI2 - rI62.*OI6;
    DI42 = phI42.*CI4 - phI24.*CI2;
    DI45 = rI45.*CI4 - rI54.*OI5;
    DI43 = phI43.*CI4 - rI34.*CI3;

    % probability of open state for InP3Rs
    OpenI = OI5 + OI6;
    % Ca increment through open InP3R channels
    DCI = krI.*OpenI.*(CeM-CcM(:,1,:));
    
    % InP3 increment from generation at plasma membrane
    DInP3Mi = kmIf.*CcM(:,3,:)./(1+kiPI*PKCv);
    % InP3 decrement, bulk degradation
    DInP3Md = kmIb*InP3M;
    
    % PKCv update
    DPKCv = kiPf*CcM(:,3,:).*(PKC0-PKCv) - kiPb*PKCv;
    
    %% RyR states
    
    rR34 = kR34.*CcMy(KRyR).^3;
    rR13 = kR13.*CcMy(KRyR).^4;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);
    
    % RyR states increments
    DR13 = phR13.*CR1 - rR31.*OR3;
    DR32 = rR32.*OR3 - rR23.*CR2;
    DR34 = phR34.*OR3 - rR43.*OR4;

    % probability of open state for RyR
    OpenR = OR3 + OR4;
    % Ca increment through open RyR channels
    DCR = krR.*OpenR.*(CeM(LRyR)-CcMy(KRyR));
        
%% Apply updates to all states %%

    %% Source/sink contributions to diffusive states (at grid centers)

    % sum of ER calcium increments
    DCeM = zeros(ngt,1);
        % from RyRs
    DCeM(LRyR) = DCeM(LRyR) - DCR;
        % from other sources (summed over all circumf. elements)
    for k=1:nc
        DCeM = DCeM + ( DCS(:,1,k) - DCI(:,1,k) - DCle(:,1,k) );
    end
    DCeM = volRce * DCeM;
    CeM = CeM + DCeM;
    
    % sum of cytosolic calcium increments
        % from RyRs
    CcMy(KRyR) = CcMy(KRyR) + DCR; % update elements with RyRs
    CcMx(:,:,1) = CcMy; % place in temp. permuted CcM array
    CcM = permute(CcMx,[1 3 2]); % permute back to standard order
        % from other ER sources
    DCcMe = -DCS + DCI + DCle;
    CcM(:,1,:) = CcM(:,1,:) + DCcMe;
        % from PM sources
    DCcMx = -DCN-DCP + DClx;
    CcM(:,3,:) = CcM(:,3,:) + DCcMx;
        % from bulk reactions
    DCcMb = DCbf-DCfb;
    CcM = CcM+DCcMb;
    
    % InP3 increments/decrements
    InP3M(:,3,:) = InP3M(:,3,:) + DInP3Mi;
    InP3M = InP3M - DInP3Md;
    
    %% diffusion processes
    
    % ER calcium diffusion (1D)
    
    %   update Ce array
    Ce = [Ce(1,1); CeM; Ce(end,1)];
    
    % axial diffusion
    DCeA = DCae*LP*Ce;
    
    % update array
    Ce = Ce + DCeA;
    % apply boundary conditions:
    Ce(1,1) = BL*Ce; % grad(Ce) (i.e. flux) = 0, distal end
    Ce(end,1) = BR*Ce; % grad(Ce) = 0, proximal end

    % Cytoplasmic calcium & InP3 diffusion (3D)
    
    %   radial diffusion
    %       temporary permuted arrays
    CcMx = permute(CcM,[2 3 1]);
    InP3Mx = permute(InP3M,[2 3 1]);
    %       perform multiplications by Laplacians
    for k=1:ngt
        DCcRx(:,:,k) = R1{k}*CcMx(:,:,k); % cyto Ca
        DInP3Rx(:,:,k) = R1{k}*InP3Mx(:,:,k); % InP3
    end
    %       permute indices back to correct order; scale by diff.coeffs.
    DCcR = DCac * permute(DCcRx,[3 1 2]); % cyto Ca
    DInP3R = DInP3 * permute(DInP3Rx,[3 1 2]); % InP3

    % circumferential diffusion
    %       temporary permuted arrays
    CcMx = permute(CcM,[3 2 1]);
    InP3Mx = permute(InP3M,[3 2 1]);
    %      perform multiplications by Laplacians
    for k=1:ngt
        for kk=1:3
            DCcCx(:,kk,k) = L1{kk,k}*CcMx(:,kk,k);
            DInP3Cx(:,kk,k) = L1{kk,k}*InP3Mx(:,kk,k);
        end
    end
    %   permute indices back to correct order; scale by diff.coeffs.
    DCcC = DCac * permute(DCcCx,[3 2 1]); % cyto Ca
    DInP3C = DInP3 * permute(DInP3Cx,[3 2 1]); % InP3
%      DCcC = zeros(ngt,3,nc);
%      DInP3C = zeros(ngt,3,nc);

    % Update arrays
    CcM = CcM + DCcR + DCcC;
    InP3M = InP3M + DInP3R + DInP3C;
    % Expand arrays to include boundaries    
    for k=1:nc
        Cc(:,:,k) = [ Cc(1,:,k); CcM(:,:,k); Cc(end,:,k) ];
        InP3(:,:,k) = [ InP3(1,:,k); InP3M(:,:,k); InP3(end,:,k) ];
    end
    
    % axial diffusion
    %       perform multiplications (need explicit loop for sparse LP)
     for k=1:nc
        DCcA(:,:,k) = LP*Cc(:,:,k); % cyto Ca
        DInP3A(:,:,k) = LP*InP3(:,:,k); % InP3
     end
     % scale by diff.coeffs.
    DCcA = DCac * DCcA;
    DInP3A = DInP3 * DInP3A;
    
    % Update Cc and InP3 arrays
    Cc = Cc + DCcA;
    InP3 = InP3 + DInP3A;
    % apply boundary conditions:
    for k=1:nc
        Cc(1,:,k) = BL*Cc(:,:,k); % grad(Cc) (i.e. flux) = 0, distal end
        Cc(end,:,k) = BR*Cc(:,:,k); % grad(Cc) = 0, prox. end
        InP3(1,:,k) = BL*InP3(:,:,k); % grad(InP3) = 0, distal end
        InP3(end,:,k) = BR*InP3(:,:,k); % grad(InP3) = 0, prox. end
        % during input period, modify  B.V.s to reflect input fluxes
        if t < tin
%             Cc(1,:,k) = Cc(1,:,k) + BVCain; % Ca pulse, distal end
%             InP3(1,:,k) = InP3(1,:,k) + BVInP3in;  % InP3 pulse, distal end
            Cc(end,:,k) = Cc(end,:,k) + BVCain; % Ca pulse, prox. end
            InP3(end,:,k) = InP3(end,:,k) + BVInP3in;  % InP3, prox. end
        end        
    end
    
    %% Calbindin state
    
    CaCalB = CaCalB + DCfb - DCbf;
    
    %% InP3R and related states
    CI1 = CI1 + DI21;
    CI2 = CI2 + DI42 - DI21 - DI26;
    CI3 = CI3 + DI43;
    CI4 = CI4 - DI45 - DI42 - DI43;
    OI5 = OI5 + DI45;
    OI6 = OI6 + DI26;
    
    PKCv = PKCv + DPKCv;
    
    %% RyR states
    CR1 = CR1 - DR13;
    CR2 = CR2 + DR32;
    OR3 = OR3 + DR13 - DR32 - DR34;
    OR4 = OR4 + DR34;
    
end

rdt = rd(ngstart:ngend);
% filename for outputs
outfile = ['OUT',mpfile];
% save data for plotting
save(outfile, 'Cct','rdt','dx','tout');
