%% This is where all output will go
TOPDIR='/Users/thompsong/Movies/Artemis';
SECONDS_PER_DAY = 86400;
addpath(fullfile(TOPDIR,'src','progressbar'))
% locate the video clips and provide their start times
origVideoPath = fullfile(TOPDIR,'VIDEO','ArtemisLaunchGT.mp4');
SECONDS_BEFORE=10;
starttime = datenum(2022,11,16,6,47,44)-SECONDS_BEFORE/SECONDS_PER_DAY;
[origVideoParent,origVideoBasename,origVideoExt] = fileparts(origVideoPath);
origVideoFramesMatDir = fullfile(TOPDIR,sprintf('origVideo%s',origVideoBasename));
figureFramesDir = fullfile(TOPDIR,'figureFrames');
mkdir(origVideoFramesMatDir);
mkdir(figureFramesDir);
% loop over each file, importing each frame and saving it as a mat-file
endtime = 0;

v=VideoReader( origVideoPath )
FPS = v.FrameRate;
s = struct('cdata',zeros(v.Height,v.Width,3,'uint8'),...
    'colormap',[]);
k = 1;
progressbar(sprintf('Capturing frames from video %s', origVideoPath));
while hasFrame(v)
    s(k).cdata = readFrame(v);
    thisframe = s(k).cdata;
    thistime = starttime + (k-1)/FPS/SECONDS_PER_DAY;
    if thistime > endtime
        endtime = thistime;
    end
    save( fullfile(origVideoFramesMatDir,sprintf('%s.mat',datestr(thistime,'HHMMSS.FFF'))), 'thisframe');
    progressbar(k/v.NumFrames) % Update progress bar
    k = k+1;
end

%
clear s c k crap thisframe thistime 

% Define the figure panel setup
close all
vw = v.Width;
vh = v.Height;
spacer = 40;
imacaspectratio=16/9;
traceHeight = 200;
figureHeight = vh + 2 * traceHeight + 4 * spacer;
figureWidth = figureHeight * imacaspectratio;

shrink_factor=4;
screen_width=1920/shrink_factor;
screen_height=1080/shrink_factor;
fh=figure('Units','pixels','Position',[0 0 screen_height*3 screen_width*3]);
screens(1)=axes(fh,'Units','pixels','Position',[0 screen_height screen_width*2 screen_height*2]); % A1-B2
screens(2)=axes(fh,'Units','pixels','Position',[0 0 screen_width screen_height]); % A3 - map?
screens(3)=axes(fh,'Units','pixels','Position',[screen_width 0 screen_width screen_height]); % B3 - seismic record section?
screens(4)=axes(fh,'Units','pixels','Position',[screen_width*2 screen_height screen_width screen_height*2]); %C1-C2 seismic spectrograms
screens(5)=axes(fh,'Units','pixels','Position',[screen_width*2 0 screen_width screen_height]); %C3 infrasound spectrograms
screens(6)=axes(fh,'Units','pixels','Position',[screen_width screen_height*2.8 screen_width*0.3 screen_height*0.1]); % text time
% plot figure files made with python in correct screens/axes
FIGDIR=fullfile(TOPDIR,'FIGURES')

imgGglEarthMap = imread(fullfile(FIGDIR,'google_earth_map3.png'));
if shrink_factor>1
    imgGglEarthMap=imresize(imgGglEarthMap,1/shrink_factor);
end
imshow(imgGglEarthMap,'Parent',screens(2));

imgSeicRecSec = imread(fullfile(FIGDIR,'artemis_record_section_seismic_ACC_Z.png'));
if shrink_factor>1
    imgSeicRecSec=imresize(imgSeicRecSec,1/shrink_factor);
end
imshow(imgSeicRecSec,'Parent',screens(3));

imgSeisSgram = imread(fullfile(FIGDIR,'specgram_seismic_ACC_scaled.png'));
imgSeisSgram=imresize(imgSeisSgram,1/shrink_factor);
imshow(imgSeisSgram,'Parent',screens(4));

% Load the seismic and infrasound data corresponding to the time window of the video files
DATADIR=fullfile(TOPDIR,'DATA')
w=waveform(fullfile(DATADIR,'acceleration.mseed'),'seed');
%plot_panels(w)
%spectrogram(w)
% create wav audio files for infrasound and seismic traces
for c=1:numel(w)
    trace_id = get(w(c),'channelinfo');
    fmmod_waveform(w(1),20+c,fullfile(TOPDIR,'AUDIO',sprintf('%s.wav',trace_id)));
end
%% plot the seismic and infrasound data
% ph1=plot(screens(4),get(w(2),'timevector'),get(w(2),'data'));
% xlim(screens(4),[starttime endtime]);
% datetick(screens(4),'x','keeplimits')
% ylabel(screens(4),'Pa')
% pylims=get(screens(4),'YLim');
% hold(screens(4));
% 
% ph2=plot(screens(5),get(w(1),'timevector'),1e-3*get(w(1),'data'));
% xlim(screens(5),[starttime endtime]);
% datetick(screens(5),'x','keeplimits')
% ylabel(screens(5),'mm/s');
% sylims=get(screens(5),'YLim');
% hold(screens(5));

%% Now loop over time from first frame to final frame
%realstart = min(starttime);
realstart = starttime;
%kmax = floor((endtime - realstart) * SECONDS_PER_DAY * FPS);
ylims3=get(screens(3),'YLim')
ylims4=get(screens(4),'YLim')
xlims3=get(screens(3),'XLim')
xlims4=get(screens(4),'XLim')
% those limits are edge of where matlab plots images. not start of axes
% within those images. add a fix factor
xlims3(1)=xlims3(1)+(xlims3(2)-xlims3(1))*0.1;
xlims3(2)=xlims3(2)-(xlims3(2)-xlims3(1))*0.1;
xlims4(1)=xlims4(1)+(xlims4(2)-xlims4(1))*0.1;
xlims4(2)=xlims4(2)-(xlims4(2)-xlims4(1))*0.1;

duration3 = 180;
duration4 = duration3;
duration1 = v.duration;
%
progressbar(sprintf('Creating ensembles frames'));
kmax = v.NumFrames;
for k=1:kmax
    progressbar(k/v.NumFrames) % Update progress bar
    relativeTime = (k-1)/FPS; % seconds passed in video
    realtime = realstart + relativeTime/SECONDS_PER_DAY;    

    xpos3 = xlims3(1) + relativeTime/duration3 * (xlims3(2)-xlims3(1));
    xpos4 = xlims4(1) + relativeTime/duration4 * (xlims4(2)-xlims4(1));
    % draw line to mark where on seismic and infrasound traces the video
    % frames being shown right now are
    
%     lh1=line(screens(3),[realtime realtime], pylims,'LineWidth',2,'Color','k');
%     lh2=line(screens(4),[realtime realtime], sylims,'LineWidth',2,'Color','k');
     lh1=line(screens(3),[xpos3 xpos3], ylims3,'LineWidth',2,'Color','k');
     lh2=line(screens(4),[xpos4 xpos4], ylims4,'LineWidth',2,'Color','k');
%     infratraveltime = 3.08/340;
%     infratime = realtime + infratraveltime/SECONDS_PER_DAY;
%     if infratime < endtime
%         lh3=line(screens(3),[infratime infratime], pylims,'LineWidth',2,'Color','r');
%         lh4=line(screens(4),[infratime infratime], sylims,'LineWidth',2,'Color','r');
%     end
%     uistack(ph1);
%     uistack(ph2);
    
    % add time
    th1=text(screens(6),0,0,datestr(realtime,'HH:MM:SS.FFF'),'FontSize',48);
    th1.Color = 'yellow';
    th1.FontSize = 48;
    % find matfiles (i.e. frames) for this time sample
    dstr = datestr(realtime,'HHMMSS.FFF');
    filepattern = fullfile(origVideoFramesMatDir,sprintf('*%s.mat',dstr));
    df=dir(filepattern);
    if numel(df)>0
        load(fullfile(df(1).folder,df(1).name));

            
        % plot the frame in the appropriate panel corresponding to posn
        % variable
        posn=1;
        if posn>0
            disp(sprintf('%s %s.%02d %d',dstr,datestr(realtime,'HHMMSS'), round(mod(realtime*SECONDS_PER_DAY,1)*30)+1,  df(1).name, posn));
            image(screens(1),thisframe)
            set(screens(1),'XTick',[],'YTick',[]);
        end


        % now all panels are plotted, save the figure window as a new jpeg
        % file
        jpgfile = fullfile(figureFramesDir,sprintf('%s.jpg',dstr));
        disp(jpgfile)
        print('-djpeg',jpgfile);
    end
    
    % delete the lines we drew to mark the time on the seismic and
    % infrasound panels
    delete(lh1)
    delete(lh2)
%     if infratime < endtime
%         delete(lh3)
%         delete(lh4)
%     end
     delete(th1)
end

%% Write video file from the JPG images - no longer need to use ImageJ which never seems to export all the images
d=dir(fullfile(figureFramesDir,'*.jpg'));
v=VideoWriter(fullfile(origVideoParent,sprintf('%s_combined.avi',origVideoBasename)),'Uncompressed AVI')
open(v)
for c=1:numel(d)
    disp(sprintf('Processing frame %d of %d',c,numel(d)))
    a=imread(fullfile(d(c).folder,d(c).name));
    writeVideo(v,a);
end
close(v);




%%
% QUICK TIME then be used to compress the 15 GB AV file to a < 1 GB .mov
% file
% iMovie can then be used to combine the .mov and .wav files, and write
% them back out to a new .mov file.

%% --------------------- FUNCTIONS FOLLOW ----------------------------- %%
function w2=fmmod_waveform(w,n,outfile)
    w=interp_waveform(w,n);
    for c=1:numel(w)
        thisw=w(c);
        fs=get(thisw,'freq')
        fc=fs*0.4;
        thisw=taper(normalize(detrend(thisw)),0.1);
        x=get(thisw,'data');
        y=fmmod(x,fc,fs,fc/4);
        w2(c)=thisw;
        w2(c)=set(w2(c),'data',y);
        %sound(y,fs)
        audiowrite(outfile,y,fs);
    end
end

function w2=interp_waveform(w,n)
    SECONDS_PER_DAY = 86400;
    for c=1:numel(w)
        x=get(w(c),'data');
        fs=get(w(c),'freq');
        t=get(w(c),'timevector');
        t2=t(1)+(1/SECONDS_PER_DAY)*(0:1/fs/n:(length(t)-1)*1/fs);
        x2=interp1(t,x,t2);
        w2(c) = w(c);
        w2(c)=set(w2(c),'data',x2);
        w2(c)=set(w2(c),'freq',fs*n);
    end
end

function w=correctUSFstations(w)
    trillium_centaur40 = 300000000.000000;
    infraBSU_centaur40 = 18.400000;
    infraBSU_centaur1 = 736.000000;
    chaparralM25_centaur40 = 160000.000000;
    chaparralM25_centaur1 = 6400000.000000;
    for c=1:numel(w)
        tr = w(c);
        Fs = get(tr,'Fs');
        net = get(tr,'network');
        sta = get(tr,'station');
        if strcmp(net, 'AM')
            continue
        end
        chan = get(tr,'channel');
        if Fs>=50
            w(c)=detrend(w(c));
            if chan(2)=='H' % a seismic velocity high-gain channel. L for low gain, N for accelerometer
                calib = trillium_centaur40;
                units = 'm/s';
            end
                                
            if chan(2)=='D' %: # infraBSU channel?
                % Assumes 0.5" infraBSU sensors at 40V p2p FS
                % But could be 1" or 5", could be 1V p2p FS or could be Chaparral M-25
                units = 'Pa';
                if (strcmp(sta,'BCHH1') & chan(3)=='1') | (strcmp(sta,'BCHH4') & chan(3)=='F')%: # Chaparral M-25
                    calib = chaparralM25_centaur40;
                elseif strcmp(sta,'BCHH') | strcmp(sta,'BCHH1') | strcmp(sta(1:3),'SAK') | strcmp(net,'NU')
                    calib = infraBSU_centaur40;
                else
                    calib = infraBSU_centaur1;
                end
            end
        
            w(c) = w(c) / calib;
            w(c) = addfield(w(c), 'calib', calib);
            w(c) = set(w(c), 'units', units);

        end
    end
end
