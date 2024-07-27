function out = videoCaCell(results_name,infile,video_name)
%%%% videoCaCell.m %%%%
% plots output of a BodyCaDynamics simulation
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
v.FrameRate = 30;

% open the video writer
open(v);

% define figure
figure(20);
ax = gca();

COLOR = colormap;

%% plot data

kdr = 5; % decimation ratio for plotting

iend = length(Ccbt); % # time steps

dph = 2*pi()/np; % azimuthal angle between processes
ns = ng; % total number of grid elements / process
ds = ld/ns; % length of grid element
es = 0.25*ds; % distance to 'protrude' process into cell (for visualization)
rb0 = rb-es; % radius at which to start depicting process (visualization)
nex = round(0.5*ns); % fraction of processes to exclude from display

dth = 2*pi()/nc; % azimuthal angle subtended by body element
nos = ncp/2; % offset to align body with processes (visualization)

lwb = 10; % line width, body elements (for visualization)
lwp = 40; % line widths, process elements (viz)

Ccmax = 4.0; % maximum expected calcium value (for colormap scaling)

for it=1:iend
    
    if ~mod(it,kdr) % decimate data to be plotted
        
        % plot body calcium
        iCcb = round( 255*Ccbt{it}/Ccmax ) + 1;
        iCcb(iCcb>256) = 256; % just in case Ccmax is exceeded
        for kr=1:nr
            r = rc(kr);
            for kc = 1:nc
                th = [ (kc-nos-1.05)*dth, (kc-nos)*dth ];
                plot( r*sin(th), r*cos(th),...
                    'Color', COLOR(iCcb(kc,kr),:),...
                    'LineWidth', lwb);
                axis square
                axis([-13 13 -13 13]);
                hold on;
            end
        end
        % plot calcium in processes
        for kp=1:np
            iCcp = mean(Ccpt{it,kp}(2:ns+1,:)')';
            iCcp = round( 255*iCcp/Ccmax ) + 1;
            iCcp(iCcp>256) = 256; % just in case Ccmax is exceeded
            ph = (kp-1)*dph;
            for ks=nex+1:ns
                sk = ns-ks;
                r = [ rb0+sk*ds, rb0+(sk+1)*ds ];
                plot( r*sin(ph), r*cos(ph),...
                    'Color', COLOR(iCcp(ks),:),...
                    'LineWidth', lwp);
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
