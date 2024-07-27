function out = videoCaCellRoot(results_name,infile,video_name)
%%%% videoCaCellRoot.m %%%%
% plots output of a BodyRootCaDynamics simulation
% contains proximal portions of input processes
% creates a video file

%% Navigate to results folder
home=pwd;
cd results
if not(isfolder(results_name))
    mkdir(results_name)
end
cd(results_name)

%% load data for plotting
load(infile);

%% create the video writer
v = VideoWriter([video_name '.avi']);
v.FrameRate = 10;

% open the video writer
open(v);

% define figure
figure(21);
ax = gca();

COLOR = colormap;

%% plot data

kdr = 5; % further decimation ratio for plotting

iend = length(Ccbt); % # time steps

% plot parameters
dph = 2*pi()/np; % azimuthal angle between processes
ns = ng+ngr; % total number of grid elements / process
ds = ld/ns; % length of grid element
es = 0.8*ds;%0.8*ds; % distance to 'protrude' process into cell (for visualization)
rb0 = rb-es; % radius at which to start depicting process (visualization)
nex = round(0.5*ns); % fraction of processes to exclude from display

dth = 2*pi()/nc; % azimuthal angle subtended by body element
nos = ncp/2; % offset to align body with processes (visualization)

lwb = 6%10; % line width, body elements (for visualization)
procradius=30;%41
lwr = procradius + procradius*0.6/(ngr*1.5)*(1:ngr); % line widths, root elements (viz)
lwp = [ procradius*ones(1,ng), lwr ];%[ 41*ones(1,ng), lwr ]; % line width, process elts. incl. root (viz)

Ccmax = 4.0; % maximum expected calcium value (for colormap scaling)

% loop over time (i.e., frames)
for it=1:iend
    
    if ~mod(it,kdr) % decimate data to be plotted
        
        % plot body calcium
        %%%   OPTION 1: Average calcium conc. over z-dimension (thickness)
            Ccb = zeros(nc,nr);
            for kt = 1:nt
                Ccb = Ccb + Ccbt{it}(:,:,kt);
            end
            Ccb = Ccb/nt;
        %%%   OPTION 2: Display calcium at a single z-value
%             kt = 5;
%             Ccb = Ccbt{it}(:,:,kt);
        % compute indices into colormap
        iCcb = round( 255*Ccb/Ccmax ) + 1;
        iCcb(iCcb>256) = 256; % just in case Ccmax is exceeded
        % plot
        for kr=1:nr
            r = rc(kr);
            for kc = 1:nc
                th = [ (kc-nos-1.05)*dth, (kc-nos)*dth ];
                plot( r*sin(th), r*cos(th),...
                    'Color', COLOR(iCcb(kc,kr),:),...
                    'LineWidth', lwb );
                axis square
                axis([-13 13 -13 13]);
                hold on;
            end
        end
        
        % plot calcium in processes
        for kp=1:np % loop over processes
            % mean over radial dimension
            iCcp = mean(Ccpt{it,kp}(2:ns+1,:)')';
            % compute indices into colormap 
            iCcp = round( 255*iCcp/Ccmax ) + 1;
            iCcp(iCcp>256) = 256;  % just in case Ccmax is exceeded
            ph = (kp-1)*dph;
            % plot
            for ks=nex+1:ns
                sk = ns-ks;
                r = [ rb0+sk*ds, rb0+(sk+1)*ds ];
                plot( r*sin(ph), r*cos(ph),...
                    'Color', COLOR(iCcp(ks),:),...
                    'LineWidth', lwp(ks) );
                axis square
                axis([-13 13 -13 13]);
                hold on;
            end
        end
        
        hold off
        
        % For making video
        frame = getframe(ax);
        frame.cdata = imresize(frame.cdata,[600, 600]);
        writeVideo(v,frame);
        disp("Another frame was saved.")
        
    end
    
end

%% close the writer object
close(v);

cd(home)
