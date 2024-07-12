% UWNoiPro = Underwater Noise Processing: main function called by 
%            graphical user interface
%
% The following file must be in the same directory as this script:
%  - evalPSDTOB[nn].m       where [nn] is current version (see below)
% 
% DESCRIPTION
% User is prompted to select one directory with all audio WAV files to be
% processed, and another one to store plots and results. A calibration
% file is required in text format with row-wise pairs of Hz-dB values.
% Processing options: sound pressure level (SPL), broadband power spectral
% density (PSD), decidecade band (DD) PSD. Processing is done in 
% intervals (snapshots) of given duration, between given start and end 
% times, or up to end of file. Snapshots may be consecutive, overlapping 
% or separated by a given time. Each snapshot is divided into segments 
% with given FFT number of points, and spectra are computed by FFT using 
% zeropad and optional overlap, to obtain an average snapshot spectrum. 
% All snapshot spectra are then averaged (arithmetic mean or median) and 
% plotted in dB units with optional high-low percentiles.
% Decidecade analysis is done in a set of consecutive bands by specifying
% start and end nominal center frequencies (ISO preferred frequencies).
% Parameters and options are passed from graphical user interface (GUI).
% Some info is output on screen during processing. Output plots are shown
% during processing for each selected option, and final ones may be
% saved to graphic files. Data may also be saved to text files.
% Calibration is passed to functions as equivalent Scale Factor in Pa.
% Decidecade bounds are computed using formula fn * 10^(m/20) with
% m = [2n-1, 2n+1] for n-th band center frequency fn (band 0 = 1 kHz). 
% DC component is subtracted from snap and an optional highpass filter
% with given cutoff frequency may be applied before processing.
% Note: Matlab R2012 or later use functions audio*(): for earlier versions
% functions wav*() are provided with different syntax. Other no longer
% recommended functions are also provided for backwards compatibility.
%
% REVISION HISTORY:
% Up to v.3.2: development versions
% v.3.3: Correct handling of first&last TOB, no illegal array index.
% v.4.1: Batch mode, param.s & opt.s declared here as static.
% v.4.2: Imported Knudsen curve overlay from 'MSPSD21.m'.
% v.4.3: Use correct file path separator for Win or Linux ('\' or '/'),
%        added option for simple/box TOB plot and max,avg in plot legend.
%        Renamed Y axis to '3rd-octave SPL' in TOB plot (was same as PSD).
%        New parameters passed to subroutines: p_ref and verb.mode.
% v.4.4: Results averaging now made on linear values, not on dB values
%        (requires broadbandspl() v.2.3, tobspl() v.2.4 or later).
% v.4.5: Changed function 'contains' to 'strcmp' due to incompatibility
%        with Matlab versions prior to R2016.
% v.4.6: Improved broadband plot with mean/median and percentiles.
%        Fixed start folder in file selection windows, and minor bugs.
% v.4.7: Revised terminology, renamed function 'broadbandspl' to 'evalpsd'
%        and 'tobspl' to 'evaltob'. Changed related filenames accordingly.
%        Added warning to increase FFT points if unable to resolve TOB.
% v.4.8: option for SPL, use mean cal.SF over full band, wave plot in Pa
%        optional highpass filter; 60 dB stopb.atten, 0.1 dB ripple, 
%        trans.band = 0.15*fc. Multiple processing types now possible, 
%        options are single boolean. Waves are always displayed first
%        in Pa using calibration, each plot is on separate figures.
% v5.02: Changed to function called by GUI developed with App Designer.
%        Parameters are passed using setup (struct type) as argument.
%        Added optional highpass filter before processing.
%        P_ref (1 microPa) no longer user defined, now static value.
%        Added command output log (diary). Return value: 0=ok, 1=error,
%        return exit status string on error. Plot titles now optional.
%        Setup optionally logged to file.
% v5.03: Plot labels, titles and legends defined globally.
% v5.04: Figures positioned near screen bottom left corner, height
%        reduced if exceeding over screen top.
% v5.05: Added dialog box for info display. Added set(currentfigure)
%        upon each plot to avoid focus stealing when clicking on another.
%        Subroutines 'evalPSd','evalTOB' changed, added version number.
%        Requires subroutines 'evalPSD' v.2.6, 'evalTOB' v.2.7.
%        Handles to current figures passed to GUI to delete them on exit.
% v5.06: Merged subroutines 'evalPSD' and 'evalTOB' to 'evalPSDBBDD' v.3.0
%        replacing one-third octave with decidecade f intervals.
%        Accepts multichannel audio files processing one channel at a time.
%
% VERSION 5.06 | 27-mar-2024 | SBU@INM

function [retVal,exStat,figHandles] = UWNoiPro506main(Setup)

VERS = '5.06'; % manually update to current version

% setup argument contains all parameters defined by user in GUI
PROCSPL     = Setup.PROCSPL;
PROCPSDB    = Setup.PROCPSDB;
PROCPSDD    = Setup.PROCPSDD;
BOXPLT      = Setup.BOXPLT;
STARTTIME   = Setup.STARTTIME;
STOPTIME    = Setup.STOPTIME;
SNAPLEN     = Setup.SNAPLEN;
SNAPOVLAP   = Setup.SNAPOVLAP;
SNAPSKIP    = Setup.SNAPSKIP;
NFFT        = Setup.NFFT;
FFTOVLAP    = Setup.FFTOVLAP;
WTYPE       = Setup.WTYPE;
AVGTYPE     = Setup.AVGTYPE;
HIPASSF     = Setup.HIPASSF;
PRCTILEL    = Setup.PRCTILEL;
PRCTILEH    = Setup.PRCTILEH;
DDIDXL      = Setup.DDIDXL;
DDIDXH      = Setup.DDIDXH;
DRAWALL     = Setup.DRAWALL;
SAVEJPG     = Setup.SAVEJPG;
SAVETXT     = Setup.SAVETXT;
XAXAUTO     = Setup.XAXAUTO;
YAXAUTO     = Setup.YAXAUTO;
XAXMIN      = Setup.XAXMIN;
XAXMAX      = Setup.XAXMAX;
YAXMIN      = Setup.YAXMIN;
YAXMAX      = Setup.YAXMAX;
LEGND       = Setup.LEGND;
FIGSIZE     = Setup.FIGSIZE;
VERB        = Setup.VERB;
PLTPAUSE    = Setup.PLTPAUSE;
KNDSDRAW    = Setup.KNDSDRAW;
KNDSSIDX    = Setup.KNDSSIDX;
% next added in v. 5.02
INPUTPATH   = Setup.INPATH; 
OUTPUTPATH  = Setup.OUTPATH;
CALFILEPATH = Setup.CALFPATH;
LOGTOFILE   = Setup.LOGFILE; % save command output to file
TITL        = Setup.TITL; % show title in all figures
LOGSETUP    = Setup.LOGSETUP; % add setup in log file
% next added in v. 5.05
KEEPFIG     = Setup.NORESIZ; % keep current figure size and positions
% end of setup argument fields

%%%%%% OLD definition of user defined parameters 
%%% Superseded by GUI input, only kept here as a reference
%%% NOTE: one-third octave definitions superseded by decidecade: however
%   center frequencies are still nominal ISO, while band limits are now
%   computed using power-of-ten formula. Extended band range [1Hz:800kHz]
%
% PROCSPL   = true;   % Evaluate SPL (true=yes, false=no)
% PROCPSD   = true;   % Evaluate PSD (true=yes, false=no)
% PROCTOB   = true;   % Evaluate TOB (true=yes, false=no)
% BOXPLT    = true;   % TOB plot type: true=boxplot, false=line plot
% STARTTIME = 0;      % Start time of analysis in s
% STOPTIME  = 30;    % Stop time of analysis in s (Inf=end of file)
% SNAPLEN   = 10;     % Snapshot duration in s
% SNAPOVLAP = 0.0;    % Snapshot overlap in percent (e.g. 50.0 = 50%)
% SNAPSKIP  = 0.0;    % Time separation between snapshots in s
% NFFT      = 2^14;   % e.g. 2^10=1024, 2^12=4096, 2^14=16384, 2^16=65536
% FFTOVLAP  = 50;     % FFT overlap in percent: typical 0 or 50
% WTYPE     = 1;      % Window: 0=none, 1=Hamming, 2=Blackman, 3=Hann
% AVGTYPE   = 1;      % Snapshot average: 1=arithmetic mean, 2=median
% HIPASSF   = 0;      % Highpass filter cutoff in Hz (0=no filter)
% PRCTILEL  = 5;      % Lower percentile (e.g. 5 = 5%, 25 = 25%)
% PRCTILEH  = 95;     % Higher percentile (e.g. 75 = 75%, 95 = 95%)
% TOBIDXL   = 4;      % First band index (see table below, as of v.5.06)
% TOBIDXH   = 31;     % Last band index  (see table below, ...)
  %   1=1Hz     11=10Hz    21=100Hz   31=1000Hz   41=10kHz    51=100kHz
  %   2=1.25 Hz 12=12.5Hz  22=125Hz   32=1250Hz   42=12.5kHz  52=125kHz
  %   ...       13=16Hz    23=160Hz   33=1600Hz   43=16kHz    ...
  %   ...       14=20Hz    24=200Hz   34=2000Hz   44=20kHz    ...
  %   ...       15=25Hz    25=250Hz   35=2500Hz   45=25kHz
  %             16=31.5Hz  26=315Hz   36=3150Hz   46=31.5kHz
  %             17=40Hz    27=400Hz   37=4000Hz   47=40kHz
  %             18=50Hz    28=500Hz   38=5000Hz   48=50kHz
  %             19=63Hz    29=630Hz   39=6300Hz   49=63kHz
  %   10=8 Hz   20=80Hz    30=800Hz   40=8000Hz   50=80kHz    60=800k
% DRAWALL   = 0;      % Draw all PSD curves (1=yes, 0=draw only low/avg/hi)
% SAVEJPG   = true;   % Save figure as *.jpg (true=yes,false=no)
% SAVETXT   = true;  % Save results as text file (true=yes,false=no)
% XAXAUTO   = true;   % PSD,TOB X axis autorange (true=yes,false=no)
% YAXAUTO   = true;   % PSD,TOB Y axis autorange (true=yes,false=no)
% XAXMIN    = 5;      % if not autorange, PSD X axis low limit in Hz
% XAXMAX    = 24000;  % if not autorange, PSD X axis high limit in Hz 
% YAXMIN    = 20;     % if not autorange, PSD,TOB Y axis low limit in dB
% YAXMAX    = 140;    % if not autorange, PSD,TOB Y axis high limit in dB 
% LEGND     = true;  % Show legend on PSD,TOB plot (true=yes,false=no)
% FIGSIZE = [1000 500]; % Size in pixels of all figures
% VERB      = 0;      % Display info (0=normal,1=detailed,2=debug)
% PLTPAUSE  = 0.5;    % Pause in s after each plot
% KNDSDRAW  = false;  % true=draw Knudsen curves, false=do not draw
% KNDSSIDX  = [0 1 2];  % sea states of Knudsen curves, 0 to 6
% PREF      = 1e-6;   % Reference pressure in Pa (standard: 1 microPa)
%%%%%% end of user defined parameters 

%%%%%% Figure labels, titles and legend - may be changed if requested
FigLab.XLABWAV  = 'Time / s';
FigLab.XLABSPL  = FigLab.XLABWAV;
FigLab.XLABPSDB = 'Frequency / Hz';
FigLab.XLABPSDD = FigLab.XLABPSDB;
FigLab.YLABWAV  = 'Sound pressure / Pa';
FigLab.YLABSPL  = 'Sound pressure level L_{p,rms} / dB re 1 \muPa';
FigLab.YLABPSDB = 'Power spectral density / dB re 1 \muPa^2/Hz';
FigLab.YLABPSDD = FigLab.YLABPSDB;
% titles are passed to sprintf() with placeholders for arguments
% note n. of placeholders for arguments: wav=3, spl=4, psd=5
FigLab.TITLWAV  = '%s: snap %d/%d';
FigLab.TITLSPL  = '%s: %.0f-%.0f s, %g s snap';
FigLab.TITLPSDB = '%s: %.0f-%.0f s, %g s snap, %d pts FFT';
FigLab.TITLPSDD = FigLab.TITLPSDB;
FigLab.LGNDPSDB = [strcat(num2str(PRCTILEL),"th percentile"); ...
                  strcat(num2str(PRCTILEH),"th percentile") ];
FigLab.LGNDPSDD = FigLab.LGNDPSDB;
% note: all strings in boxplot legend must be equal length
FigLab.LGNDBOX = ['Box edges = 25th-75th percentiles'; ...
                  'Box center line = median         '; ...
                  'Dashed vertical line = min-max   '; ...
                  'Cross = outlier                  ' ];

%%%%%% static global parameters - do not change
PREF      = 1e-6;   % Reference pressure in Pa (standard: 1 microPa)
ISOD = [1 1.25 1.6 2 2.5 3.15 4 5 6.3 8]; % ISO preferred f in a decade
ISOF = [ISOD ISOD*10 ISOD*100 ISOD*1e3 ISOD*1e4 ISOD*1e5]; % 1Hz - 800kHz
% Knudsen's reference noise levels for sea states 0 to 6 (Ref: Lurton)
KNDLEVS = [44.5  55  61.5  64.5  66.5  68.5  70]; % @1kHz, SS0 to SS6
SSLABELS = ['SS0';'SS1';'SS2';'SS3';'SS4';'SS5';'SS6']; % Knudsen labels
KNDDEC = (-17);   % NL = KNDLEVS + KNDDEC*log10(f in kHz) for f > KNDCUT
KNDCUT = 1e3;     % level assumed constant up to this frequency in Hz
TXTXPOS = 1.5*XAXMIN;   % X position for knudsen curve labels
TXTYOFFS = 1.5;         % label Y offset in dB with respect to line
FIGWAV = 1; FIGSPL = 2; FIGPSDB = 3; FIGPSDD = 4; % figure handles
%FIGINFO = 9; % additional figure for info display during processing
DISPINFO = true; % true=show info figure, false=do not show
FIGXYPOS = [60 120]; % initial X Y fig.pos., pix from bottom left corner
FIGXYOFFS = [50 -30]; % incremental X Y offset for each figure in pix
FILTBOXPOS = [0.3 0.3 0.4 0.1]; % "filtering..." box msg. pos.: X Y W H
FIGDRV = '-djpeg'; % graphics driver for figure output (e.g. jpeg)
FIGRES = '-r300'; % output figure resolution in dpi
TXTDECS = 4; % n. of decimals in output text data (default = 13)
LOGFILENAME = 'UWNoiPro_log.txt'; % log file to save command window output
%%%%%% end of static global parameters 

if ispc  % file path separator: '\' for Windows, '/' for Linux
    PathSeparator = '\';
else
    PathSeparator = '/';
end
matlabVers = version('-release');
versYear = str2num( matlabVers(1:4) ); % skip trailing 'a' or 'b'

retVal = 0; % clear return value to 0 (no error) at beginning
exStat = ''; % on error, return error description string
figHandles = ''; % do not return figure handles on error

timeOfStart = datetime('now');          
if LOGTOFILE
    diary([OUTPUTPATH,PathSeparator,LOGFILENAME])
    disp(char(61*ones(1,80))) % separator line: ascii 61 = '='
    fprintf('Start: %s\n',timeOfStart)
else
    diary off
end

% display header in command window and logfile if enabled
hdr = ['| UWNoiPro = Underwater Noise Processing v. ',VERS,' |'];
hdrb = char([111 45*ones(1,length(hdr)-2) 111]); % ascii 45='-', 111='o' 
disp([hdrb; hdr; hdrb]);

%%%%% select folders input, output, and calib.file
InputFolder = INPUTPATH;
OutputFolder = OUTPUTPATH;
wavPath = [InputFolder PathSeparator];  
outPath = [OutputFolder PathSeparator];
wFiles = dir(InputFolder);  

validWFiles = 0;
for idx = 3:length(wFiles)   % count all valid WAV files skipping current
                             % and parent directory which come first
    len = length(wFiles(idx).name);
    if len > 3
        ext = lower(wFiles(idx).name(len-3:len)); % .wav for valid file
    else
        ext = [];   % not a valid filename with 3-char extension
    end
    if wFiles(idx).isdir == 0 && strcmp(ext,'.wav')
        validWFiles = validWFiles+1;
    end
end
fprintf('Input folder %s:\n  found %d WAV files out of %d total\n',...
            InputFolder,validWFiles,length(wFiles)-2);
fprintf('Output folder: %s\n', OutputFolder);
        
%%%%% load calibration file
CalFilePath = CALFILEPATH;
if CalFilePath == 0 % user canceled
    exStat = 'Error: must select a valid calibration file.';
    disp(exStat)
    retVal = 1;
    msgbox(exStat,'Error','error','modal');
    return
else
    if ~isfile(CalFilePath) % initial condition, no file
        exStat = 'Error: missing calibration file.';
        disp(exStat)
        retVal = 1;
        msgbox(exStat,'Error','error','modal');
        return
    end
end
[~,CalFileName,CalFileExt] = fileparts(CalFilePath); % folder is unused
calTab = readtable(CalFilePath); % blank lines give NaN
calMat = calTab{:,:}; % convert table into matrix for further processing
calMat = rmmissing(calMat); % remove NaN if blank lines in table
calPts = size(calMat,1);
if calPts < 3
    exStat = 'Error: at least 3 points required in calibration file.';
    disp(exStat)
    retVal = 1;
    msgbox(exStat,'Error','error','modal');
    return
end

if max(calMat(:,2)) > 0 % assume file is in SF Pa units, not dB
    msg = 'Calibration file seems in SF(Pa) units, converting to dB';
    disp(msg)
    msgbox(msg,'Warning','warn','modal');
    % convert SF values to dB
    calMat(:,2) = 20*log10( 1/calMat(:,2) ) - 120;
end
fprintf('Calibration file %s:\n  %d points spanning (%g-%g) Hz\n',...
    CalFilePath, calPts, calMat(1,1), calMat(calPts,1) );

refdb = 20 * log10(PREF); % all levels are relative to reference p

%%%%%%%%%%%%%%%%%%%%% main loop over audio files %%%%%%%%%%%%%%%%%%%%%%%%

for wFileN = 3:length(wFiles)

    len = length(wFiles(wFileN).name);
    if len > 3
        ext = lower(wFiles(wFileN).name(len-3:len)); % .wav for valid file
    else
        ext = [];   % no filename extension
    end

    if wFiles(wFileN).isdir > 0 || ~strcmp(ext,'.wav')
            fprintf('\n"%s": not a WAV file\n',wFiles(wFileN).name)
        continue    % skip to next audio file
    else
        % valid file, proceed with processing
        WavFileName = wFiles(wFileN).name;
    end
    fprintf('\nInput file %d of %d: "%s"\n',...
                        wFileN-2, length(wFiles)-2, WavFileName);

    % initialize figures, later on to be called by set(0,'currentfigure')
    figWav=figure(FIGWAV); set(gcf,'Name','Wave');
    if PROCSPL
        figSpl=figure(FIGSPL); set(gcf,'Name','SPL');
    end
    if PROCPSDB
        figPsdB=figure(FIGPSDB); set(gcf,'Name','Broadband PSD');
    end
    if PROCPSDD
        figPsdD=figure(FIGPSDD); set(gcf,'Name','Decidecade PSD');
    end

    WavFile = [wavPath WavFileName]; % not affected by trailing '(chN)'
                                     % if multichannel

    if versYear < 2012
    %    function fileparts() returns different values prior to R2012
        [F.path,F.name,F.ext,F.vers] = fileparts(WavFile);
    else
        [F.path,F.name,F.ext] = fileparts(WavFile);
    end
    FileNameNoExt = F.name; % output files share same filename

    if versYear < 2012
    % audio file info: limited support prior to R2012a, need to use
    % function 'wavinfo' which also returns additional unused info
        [~, info2] = wavfinfo(WavFile); % info strings, 1st unused
        % need to parse return string to get total samples
        WavFileInfo.TotalSamples = sscanf(info2,'%*s %*s %*s %*s %d');
        [~,p2,~] = wavread(WavFile,1); % read 1 smp to get srate
        WavFileInfo.SampleRate = p2;
        WavFileInfo.Duration = WavFileInfo.TotalSamples / p2;
    else
        WavFileInfo = audioinfo(WavFile); % contains all required info
    end
    smpRate = WavFileInfo.SampleRate;
    nSmp    = WavFileInfo.TotalSamples;
    fileDur = WavFileInfo.Duration; % in seconds
    fileDurStr = duration(0,0,fileDur,'Format','hh:mm:ss.SSS');
    nChans  = WavFileInfo.NumChannels;
    if nSmp > 0
        fprintf(' Rate: %g kHz, %.3f Msamples, duration: %s (%.3f s)\n',...
            smpRate/1e3, nSmp/1e6, fileDurStr, fileDur);
        if not(isinf(STOPTIME)) && STOPTIME > fileDur
            exStat = '* Warning: stop time exceeds file duration.';
            disp(exStat)
            retVal = 1;
            validWFiles = validWFiles - 1; % files actually processed
            continue % does not exit program, skip to next file
        end
    else
        exStat = '* Warning: acoustic file seems to be void.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % does not exit program, skip to next file
    end

    WavFileNameOrig = WavFileName; % to append (chN) if multichannel
    FileNameNoExtOrig = FileNameNoExt;
    
    %%%%%%%%%%%%%%%%% loop over channels if more than 1 %%%%%%%%%%%%%%%

    for curChan = 1:nChans 
        % loop index also used to avoid repeating global info display

    if nChans > 1 % stereo or multichannel, append (chN) to filename
        WavFileName   = [ WavFileNameOrig   sprintf('(ch%d)',curChan)];
        FileNameNoExt = [ FileNameNoExtOrig sprintf('(ch%d)',curChan)];
    end

    if isinf(STOPTIME)
        stopTime = fileDur;
    else
        stopTime  = STOPTIME;
    end
    snapLen   = SNAPLEN;
    snapOvlap = SNAPOVLAP;
    snapSkip  = SNAPSKIP;

    if snapLen > fileDur
        exStat = '* Warning: snapshot exceeds file duration.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % does not exit program, skip to next file
    end

    % manage snapshot overlap and skip time
    nSnapOvlap = (100/(100-snapOvlap)); % snap.s starting within snapLen
    % nSnap is total no. of snapshots accounting for possible overlap
    nSnapExact = nSnapOvlap*((stopTime-STARTTIME)/snapLen - 1);
    nSnap = floor(nSnapExact) + 1;
    snapStep = snapLen / nSnapOvlap; % time step between snapshots
    % override if skip > 0 (get rid of overlap, should be 0 anyway)
    if snapSkip > 0
        nSnap = floor((stopTime-STARTTIME)/(snapLen+snapSkip));
        snapStep = snapLen + snapSkip;
    end

    % build matrix of [start,stop] samples for each snap
    snapSmp  = zeros(nSnap,2); % 1st col=start, 2nd col=stop
    for idx = 1:nSnap % set start & stop smp for each snap
        snapSmp(idx,1) = round((STARTTIME+(idx-1)*snapStep)*smpRate+1);
        snapSmp(idx,2) = round(snapSmp(idx,1) + snapLen*smpRate - 1);
    end
    snapTime = snapSmp/smpRate; % start & stop time for each snap

    % define snapshot timebase
    Snap.T = 0 : 1/smpRate : snapLen-(1/smpRate); % struct: Snap

    if curChan == 1 % do not repeat info display for next channels
        fprintf(' Timebase: [%.1f-%.1f] s, %d snapshots, %.1f s each', ...
            STARTTIME,stopTime,nSnap,snapLen);
        if SNAPOVLAP > 0
            fprintf(', %d%% overlap', SNAPOVLAP);
        end
        if SNAPSKIP > 0
            fprintf(', separated by %.1f s', SNAPSKIP);
        end
        fprintf('\n');
    end

    % manage FFT overlap
    nfftOvlap = (100/(100-FFTOVLAP)); % FFTs starting within each nFFT
    % n. of FFTs in snapshot, accounting for possible overlap:
    snapFFTnexact = nfftOvlap*(snapLen*smpRate/NFFT - 1) + 1;
    snapFFTn = ceil(snapFFTnexact); % round up for zeropad
    offsFFT = NFFT / nfftOvlap; % incremental offset of each FFT
    % total samples spanned by FFT, accounting for possible overlap:
    spanFFT = NFFT + (snapFFTn-1)*offsFFT; 

    if snapFFTnexact-fix(snapFFTnexact) > 0 % noninteger FFTn: zero pad
        zeroPad = spanFFT - snapLen*smpRate;
    else   % exact n. of FFTs in snapshot: no zero pad
        zeroPad = 0;  
    end

    % broadband f vector ranges from f=0 to fmax=(smp.rate)/2
    fVect = linspace(0, smpRate/2, NFFT/2 + 1);

    % Interpolate calibration data (in dB) in FFT points,
    % use 'nearest' method to keep flat response outside range
    % NOTE: to be done for each file as it depends on sample rate
    CaldB = interp1(calMat(:,1),calMat(:,2),fVect,'nearest','extrap');
    CalSF  = 1./(10.^((CaldB+120)/20)); % Scale Factor in Pa
    meanCalSF = mean(CalSF);  % for waveform and SPL

    wType = WTYPE;
    switch wType
        case 0
            FFTW = ones(NFFT,1); % rectangular (no window)
        case 1
            FFTW = hamming(NFFT);
        case 2
            FFTW = blackman(NFFT);
        case 3
            FFTW = hann(NFFT);
        otherwise
            FFTW = ones(NFFT,1);
            if VERB > 0
                disp('* Invalid window, using rectangular');
            end
    end
    avgType = AVGTYPE;  % arithmetic mean or median

    % set of parameters to be passed on to external functions
    Param.sr = smpRate;
    Param.nf = NFFT;    % FFT n. of points 
    Param.av = avgType;
    Param.fn = snapFFTnexact; % n. of FFTs in each snap
    Param.of = offsFFT;  % offset between FFTs
    Param.zp = zeroPad;  % zero padding for last FFT
    Param.re = PREF;    % reference pressure
    Param.vb = VERB;    % verbose mode
    Param.pp = PLTPAUSE; % pause after each plot
    Param.procpsdb = PROCPSDB; % boolean
    Param.procpsdd = PROCPSDD; % boolean
    Param.figwav = figWav; % wav figure handle
    if PROCPSDB
        Param.figpsdb = figPsdB; % broadband psd plot figure handle
    end
    if PROCPSDD
        Param.figpsdd = figPsdD; % decidecade psd plot figure handle
    end

    if curChan == 1 % do not repeat info display for next channels
        if PROCPSDB || PROCPSDD
            fprintf(' FFT: %d points, %d%% overlap', Param.nf,FFTOVLAP);
            fprintf(', resolution = %.3f Hz\n', fVect(2));
            fprintf(' Snapshot contains %g FFT segments\n', Param.fn);
        end
    end

    %%%%%%% Decidecade (DD) band setup: center f and limits
    if PROCPSDD
        DDmat = []; % initialize DD index matrix
        for idx = DDIDXL:DDIDXH
            % Decidecade center f are nominal ISO preferred, while
            % band limits are computed using power-of-ten formula.
            % Index 1 is for 1-Hz band, need to shift by 1 to get 10^0.
            DDfl = 10^( (2*(idx-1)-1) / 20 );
            DDfh = 10^( (2*(idx-1)+1) / 20 );
            fVidxl = find(fVect>=DDfl, 1, 'first'); % fvect low idx
            fVidxh = find(fVect<=DDfh, 1, 'last');  % fvect hi idx
            if DDfl > smpRate/2 || DDfh > smpRate/2
                exStat = 'Error: topmost decidecade too high.';
                disp(exStat)
                retVal = 1;
                msgbox(exStat,'Error','error','modal');
                return
            end
            if VERB > 0
                fprintf('DD(%.1f)=[%.3f-%.3f], [%d-%d]=[%.3f-%.3f]\n', ...
                    ISOF(idx),DDfl,DDfh, ...
                    fVidxl, fVidxh, fVect(fVidxl),fVect(fVidxh) );
            end           
            if fVect(2) > DDfh - DDfl % f resolution too coarse
                exStat='Error: increase FFT points to resolve DD band.';
                disp(exStat)
                retVal = 1;
                msgbox(exStat,'Error','error','modal');
                return
            end
            %%% DD matrix: ISO [f0, flow, fhi]; computed [flow, fhi];
            %   vector indexes [low, hi]
            DDmat = [DDmat; ISOF(idx) DDfl DDfh ...
                fVect(fVidxl) fVect(fVidxh) fVidxl fVidxh];
        end
        DDf0 = DDmat(:,1); % vector of ISO preferred center freq.s
        % interpolate calibration SF values in ISO center frequencies
        DDCalSF = interp1(fVect,CalSF,DDf0,'linear');
    else
        % pass void arguments to subroutine if DD not selected
        DDmat = []; DDCalSF = []; 
    end %%%%%%%%%% Decidecade band setup
    
    % initialize each results matrix
    outMatSpl = []; outMatPsdB = []; outMatPsdD = [];
    
    % reduce figure height for small screens
    scrSiz = get(0,'ScreenSize');
    if FIGSIZE(2)+FIGXYPOS(2) > scrSiz(4) % figure top goes out of screen
        FIGSIZE(2) = scrSiz(4)-FIGXYPOS(2);
    end

    % clear and reset position of all figures if requested
    clf(figWav);
    if PROCSPL
        clf(figSpl)
    end
    if PROCPSDB
        clf(figPsdB)
    end
    if PROCPSDD
        clf(figPsdD)
    end
    if ~KEEPFIG
        set(0,'CurrentFigure',figWav);
        figWav.OuterPosition=[FIGXYPOS FIGSIZE];
        if PROCSPL
            set(0,'CurrentFigure',figSpl);
            figSpl.OuterPosition=[FIGXYPOS+FIGXYOFFS FIGSIZE];
        end
        if PROCPSDB
            set(0,'CurrentFigure',figPsdB);
            figPsdB.OuterPosition=[FIGXYPOS+2*FIGXYOFFS FIGSIZE];
        end
        if PROCPSDD
            set(0,'CurrentFigure',figPsdD);              
            figPsdD.OuterPosition=[FIGXYPOS+3*FIGXYOFFS FIGSIZE];
        end
    end
    if DISPINFO % set info figure contents and properties
        infoTx1=sprintf(['WAV: %s\n    sr = %.1f kHz,    ' ...
            'len = %.1f s\n\n'], WavFileName, smpRate/1e3, fileDur);
        infoTx2=sprintf(['CAL: %s\n    %d pts, (%g-%g) Hz,    ' ...
            'mean SF = %g\n\n'], [CalFileName CalFileExt], ...
            calPts, calMat(1,1), calMat(calPts,1), meanCalSF );
        infoTx3=sprintf('FFT res. = %.2f Hz\n\n', fVect(2) );
        if nSnap > 1
            infoTx4=sprintf('%d snaps,    %.2f FFT segments each\n', ...
                nSnap, Param.fn);
        else
            infoTx4=sprintf('1 snap,    %.2f FFT segments\n', Param.fn);
        end
        figInfo = msgbox([infoTx1 infoTx2 infoTx3 infoTx4], ...
            'Info','replace'); % redraws on same box for each file
        figInfoPos = get(figInfo,'OuterPosition');
        set(figInfo,'OuterPosition',[0 0 figInfoPos(3) figInfoPos(4)]);
    end

    if nChans > 1 % stereo or multichannel
        fprintf(' Processing channel %d of %d:\n',curChan,nChans);
    end

    %%%%%%%%%%%%%%%%%%%%%%% snapshot loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for snapN = 1 : nSnap

        Param.sn = snapN;  % passed on to subroutines

        smpToRead = snapSmp(snapN,:); % 2-element vect [start,stop]
        fprintf(' Snap %d/%d t=[%.1f-%.1f] s=[%d-%d] ',...
            snapN,nSnap,smpToRead/smpRate,smpToRead);
        
        % read snapshot from file: use different function prior to R2012
        %%%%% WARNING: reading multichannel with 'wavread()' not checked
        if versYear < 2012
            [AllChanY,smpRate,~]=wavread(WavFile,smpToRead);
            % 3rd return value above (nBits) is unused
        else
            [AllChanY,smpRate] = audioread(WavFile,smpToRead);
        end

        if nChans > 1
            OrigY = AllChanY(:,curChan); % select current channel
        else
            OrigY = AllChanY; % mono
        end

        Snap.Y = OrigY - mean(OrigY); % remove DC component
        snapStart = smpToRead(1)/smpRate; % update start time

        %%%%% plot snap waveform first anyway
        set(0,'CurrentFigure',figWav);
        plot(Snap.T + snapStart, Snap.Y * meanCalSF)
        grid on
        xlabel(FigLab.XLABWAV)
        ylabel(FigLab.YLABWAV)
        if TITL
            title(sprintf(FigLab.TITLWAV,WavFileName,snapN,nSnap), ...
              'Interpreter','none') % otherwise underscores in filename
                                    % are interpreted as subscript
        else
            title('')
        end
        if HIPASSF > 0
            % draw annotation, filter, replot filtered wave
            ann = annotation('textbox',FILTBOXPOS, ...
                'string','Filtering ...');
            ann.HorizontalAlignment = 'center';
            ann.VerticalAlignment = 'middle';
            ann.BackgroundColor = [1 1 1];
            fprintf('filtering... '); % may take some time to filter
            pause(PLTPAUSE)
            Snap.Y = highpass(OrigY, HIPASSF, smpRate); 
            % default atten. 60dB, ripple 0.1dB, trans.width = 15%*f_cut
            clf; 
            plot(Snap.T + snapStart, Snap.Y * meanCalSF)
            grid on % re-apply grid, labels and title as above
            xlabel(FigLab.XLABWAV)
            ylabel(FigLab.YLABWAV)
            if TITL
                title(sprintf([FigLab.TITLWAV ' filtered'], ...
                    WavFileName,snapN,nSnap), 'Interpreter','none')
            else
                title('')
            end
        end
        pause(PLTPAUSE) 
        
        if PROCSPL 
            %%%%%%%%%%% SPL: square of snap RMS by mean SF
            snapSpl = (rms(Snap.Y) * meanCalSF)^2;
            outMatSpl = [outMatSpl snapSpl]; % it is actually a vector
            outdbMatSpl  = 10*log10(outMatSpl)  - refdb;
            % update SPL RMS plot with current snap data 
            set(0,'CurrentFigure',figSpl);
            % points are drawn along x at snapshot midpoints
            snapMidTimes = (snapTime(1:snapN,2)+snapTime(1:snapN,1)) / 2; 
            plot(snapMidTimes,outdbMatSpl,':ob','MarkerFaceColor','b')
            grid on
            xlabel(FigLab.XLABSPL)
            ylabel(FigLab.YLABSPL)
            ax = gca;  % get current axes handle
            ax.XScale = 'linear';
            xlim([STARTTIME stopTime]) % set time axis to full interval
            if YAXAUTO
                ylim auto % autorange: override Y limits
            else
                ylim([YAXMIN YAXMAX])
            end
            if TITL
                title(sprintf(FigLab.TITLSPL, ...
                    WavFileName,STARTTIME,nSnap*snapLen,snapLen), ...
                    'Interpreter','none') 
            else
                title('')
            end
        end
        pause(PLTPAUSE)

        %%%%%%%%%%%% PSD subroutine call, processing in BB and/or DD

        if PROCPSDB || PROCPSDD 
            [curPsdB,curPsdD] = ...
                    evalPSDBBDD30(Snap,FFTW,CalSF,DDCalSF,DDmat,Param);

            outMatPsdB = [outMatPsdB curPsdB];
            outMatPsdD = [outMatPsdD curPsdD];
        end
        % pause already included at end of subroutine

        fprintf('\n'); % newline in command window when done

    end % snapshot loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% BB, DD averaging (done on linear values, not dB)
    switch avgType 
        case 1
            avgStr = 'Arithmetic mean';
            outAvgPsdB = mean(outMatPsdB,2);
            outAvgPsdD = mean(outMatPsdD,2);
        case 2
            avgStr = 'Median';
            outAvgPsdB = median(outMatPsdB,2);
            outAvgPsdD = median(outMatPsdD,2);
        otherwise % default: should not get here anyway
            avgStr = 'Arithmetic mean';
            outAvgPsdB = mean(outMatPsdB,2);
            outAvgPsdD = mean(outMatPsdD,2);
    end

    %%%%% percentiles for psd
    outPctHPsdB = prctile(outMatPsdB,PRCTILEH,2);
    outPctLPsdB = prctile(outMatPsdB,PRCTILEL,2);
    outPctHPsdD = prctile(outMatPsdD,PRCTILEH,2);
    outPctLPsdD = prctile(outMatPsdD,PRCTILEL,2);

    %%%%% convert results to dB re 1 microPa
    outdbMatPsdB   = 10*log10(outMatPsdB)  - refdb; % full matrix 
    outdbAvgPsdB   = 10*log10(outAvgPsdB)  - refdb; % average vector
    outdbPctHPsdB  = 10*log10(outPctHPsdB) - refdb;
    outdbPctLPsdB  = 10*log10(outPctLPsdB) - refdb;
    outdbMatPsdD   = 10*log10(outMatPsdD)  - refdb; % full matrix 
    outdbAvgPsdD   = 10*log10(outAvgPsdD)  - refdb; % average vector
    outdbPctHPsdD  = 10*log10(outPctHPsdD) - refdb;
    outdbPctLPsdD  = 10*log10(outPctLPsdD) - refdb;

    %%%%%%%%%%%%%%%%%%%%%% plot results according to processing type
    if PROCSPL
        % no further plot needed: figure was updated for each snap
    end
    
    if PROCPSDB
        set(0,'CurrentFigure',figPsdB);
        clf % do not redraw on last snapshot plot
        hold on
        if DRAWALL > 0
            for snapN = 1 : nSnap % draw all snapshots
               semilogx( fVect, outdbMatPsdB(:,snapN), 'c' )
            end
        end
        % return values are used to reference legend items
        psdL = semilogx( fVect, outdbPctLPsdB, 'g' ); % low percentile
        psdH = semilogx( fVect, outdbPctHPsdB, 'r' ); % hi percentile
        psdA = semilogx( fVect, outdbAvgPsdB, 'b' ); % avg: on top
        grid on
        hold off
    end
    
    if PROCPSDD
        % line plot: selected average type
        % boxplot: median, fixed 25%-75% percentiles, min-max, outliers
        % transpose matrix, series must be column-wise
        set(0,'CurrentFigure',figPsdD);
        if BOXPLT 
            if nSnap == 1 % no statistics for percentiles, use line plot
                fprintf('\n Single snapshot: not enough for boxplot\n');
                semilogx(DDf0,outdbAvgPsdD,'--ob', ...
                    'MarkerFaceColor','b') % plot avg line with points
            else % plot both avg line and boxplot with percentiles
                plot(1:length(DDf0),outdbAvgPsdD,'--b')
                hold on
                boxplot(outdbMatPsdD','Labels',num2str(DDf0))
                hold off
            end
        else % plot three lines for low-hi percentiles and average
            ddL = semilogx(DDf0,outdbPctLPsdD,'og');
            hold on
            ddH = semilogx(DDf0,outdbPctHPsdD,'or');
            ddA = semilogx(DDf0,outdbAvgPsdD,'ob','MarkerFaceColor','b');
            hold off
        end
    end

    %%%%% overlay Knudsen curves for N sea states - only for psd
    if PROCPSDB && KNDSDRAW
        set(0,'CurrentFigure',figPsdB);
        hold on
        seaStates = SSLABELS; 
        kndsIdx = KNDSSIDX + 1; % SS range from 0, index from 1
        for idx = 1:length(kndsIdx)
            knudsNL1K = KNDLEVS(kndsIdx(idx));
            knudsNL = knudsNL1K + KNDDEC*log10(fVect/1e3); 
            txtYpos = knudsNL1K + KNDDEC*log10(fVect(1)/1e3); 
            for idx2 = 1:length(fVect)
                if fVect(idx2) < KNDCUT
                    knudsNL(idx2) = knudsNL1K;
                    txtYpos = knudsNL1K;
                end
            end
            semilogx( fVect, knudsNL, '--m' )
            text(TXTXPOS,txtYpos+TXTYOFFS,seaStates(kndsIdx(idx),:))
        end
        hold off
    end

    %%%%%%%%% apply type-specific title, scale, labels, opt. legend.
    % Note: interpreter set to 'none' in title otherwise underscores
    %       in filename are interpreted as subscript.

    if PROCSPL
        % plot already complete
    end

    if PROCPSDB
        set(0,'CurrentFigure',figPsdB);
        grid on
        xlabel(FigLab.XLABPSDB)
        ylabel(FigLab.YLABPSDB)
        ax = gca;  % get current axes handle
        ax.XScale = 'log';
        if XAXAUTO
            xlim auto % autorange: override X limits
        else
            xlim([XAXMIN XAXMAX])
        end
        if YAXAUTO
            ylim auto % autorange: override Y limits
        else
            ylim([YAXMIN YAXMAX])
        end
        if TITL
            title(sprintf(FigLab.TITLPSDB,...
                WavFileName,STARTTIME,nSnap*snapLen,snapLen,NFFT), ...
                'Interpreter','none')
        else
            title('')
        end
        if LEGND
            legend( [psdH psdA psdL], ...
                FigLab.LGNDPSDB(2), avgStr, FigLab.LGNDPSDB(1) )
        end
    end

    if PROCPSDD
        set(0,'CurrentFigure',figPsdD);
        grid on
        xlabel(FigLab.XLABPSDD)
        ylabel(FigLab.YLABPSDD)
        ax = gca;  % get current axes handle
        if BOXPLT
            if nSnap == 1
                ax.XScale = 'log'; % no stats for percentiles
            else
                ax.XScale = 'linear'; % boxplot decidecade series
            end
        else
            ax.XScale = 'log'; % decidecade line plot
        end
        ylim([YAXMIN YAXMAX]) % only resize Y axis, X set by DD bands
        if YAXAUTO
            ylim auto % autorange: override Y limits
        end
        if TITL
            title(sprintf(FigLab.TITLPSDD,...
                WavFileName,STARTTIME,nSnap*snapLen,snapLen,NFFT), ...
                'Interpreter','none') 
        else
            title('')
        end
        if LEGND
            if BOXPLT && nSnap > 1 % No boxplot output for single snap.
                % Use trick to split legend in multiple rows to
                % align linetype mark with string of average type;
                % rows must be equal length, pad avgStr with spaces
                spPad = length(FigLab.LGNDBOX)-strlength(avgStr);
                legend([ FigLab.LGNDBOX(1,:); ...
                         FigLab.LGNDBOX(2,:); ...
                         [avgStr char(32*ones(1,spPad))]; ...
                         FigLab.LGNDBOX(3,:); ...
                         FigLab.LGNDBOX(4,:) ]) 
            end
            if ~BOXPLT % decidecade line plot, also for single snap
                legend( [ddH ddA ddL], ...
                FigLab.LGNDPSDD(2), avgStr, FigLab.LGNDPSDD(1) )

            end
        end
    end

    % Add date_time in output filenames (jpeg&text) to make them unique
    dtStr = string(datetime('now'), '_yyMMdd_HHmm');
    
    %%%%% Save figures as *.jpg file with same size as on screen
    fullPath = [outPath FileNameNoExt]; % File path and name without ext.
    if SAVEJPG
        disp([' Saving figures to ' outPath ' :'])
        if PROCSPL
            set(FIGSPL,'PaperPositionMode','auto')
            jpegFile = strcat(fullPath, dtStr, '_spl.jpg');
            print(['-f' num2str(FIGSPL)],FIGDRV,FIGRES,jpegFile)
            [~,fName,fExt] = fileparts(jpegFile);
            disp(strcat(' - ', fName, fExt))
        end
        if PROCPSDB
            set(FIGPSDB,'PaperPositionMode','auto')
            jpegFile = strcat(fullPath, dtStr, '_psd.jpg');
            print(['-f' num2str(FIGPSDB)],FIGDRV,FIGRES,jpegFile)
            [~,fName,fExt] = fileparts(jpegFile);
            disp(strcat(' - ', fName, fExt))
        end
        if PROCPSDD
            set(FIGPSDD,'PaperPositionMode','auto')
            if BOXPLT
                jpegFile = strcat(fullPath, dtStr, '_ddbox.jpg');
                print(['-f' num2str(FIGPSDD)],FIGDRV,FIGRES,jpegFile)
                [~,fName,fExt] = fileparts(jpegFile);
                disp(strcat(' - ', fName, fExt))
            else
                jpegFile = strcat(fullPath, dtStr, '_dd.jpg');
                print(['-f' num2str(FIGPSDD)],FIGDRV,FIGRES,jpegFile)
                [~,fName,fExt] = fileparts(jpegFile);
                disp(strcat(' - ', fName, fExt))
            end
        end
    end

    %%%%%% Save results as text files.
    % Transpose results vectors to build output matrix,
    % use rounding to reduce n. of decimals from default = 13.
    % Note: function dlmwrite() not recommended after R2019,
    %       use function writematrix() instead
    if SAVETXT
        disp([' Saving text files to ' outPath ' :'])
        if PROCSPL
            outTxtMat = [snapMidTimes,outdbMatSpl'];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_spl.txt');
            if versYear < 2019 
                dlmwrite(outTxtFile,outTxtMat,'\t') % not recommended
            else
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(strcat(' - ', fName, fExt))
        end

        if PROCPSDB
            outTxtMat = [fVect',outdbAvgPsdB,outdbPctLPsdB,outdbPctHPsdB];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_psd.txt');
            if versYear < 2019 
                dlmwrite(outTxtFile,outTxtMat,'\t') % not recommended
            else
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(strcat(' - ', fName, fExt))
        end

        if PROCPSDD
            outTxtMat = [DDf0,outdbAvgPsdD,outdbPctLPsdD,outdbPctHPsdD];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_dd.txt');
            if versYear < 2019 
                dlmwrite(outTxtFile,outTxtMat,'\t') % not recommended
            else
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(strcat(' - ', fName, fExt))
        end
    end
    pause(PLTPAUSE * 4) % processing end, longer pause

    end % channels loop

end %%%%%%%%%% datafile loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pass existing figure handles to GUI so they can be closed on exit
if exist('figWav','var') > 0
    figHandles = figWav; % Should exist unless an error occurs
end
if exist('figInfo','var') > 0
    figHandles = [figHandles figInfo];
end
if exist('figSpl','var') > 0
    figHandles = [figHandles figSpl];
end
if exist('figPsdB','var') > 0
    figHandles = [figHandles figPsdB];
end
if exist('figPsdD','var') > 0
    figHandles = [figHandles figPsdD];
end

timeOfEnd = datetime('now');
elapsTime = sprintf(' in %s', timeOfEnd - timeOfStart);
switch validWFiles
    case 1
        exStat = 'Normal end - Processed 1 file';
    otherwise
        exStat = sprintf('Normal end - Processed %d files',validWFiles);
end
exStat = [exStat elapsTime]; % Add elapsed processing time

if LOGTOFILE
    if LOGSETUP % Display setup so it can also be logged to file
        fprintf('\nSetup:\n');
        disp(Setup)
    end
    fprintf('\nDone: %s\n\n',datetime('now'))
    diary off % Stop logging to file before exit, if enabled
else
    disp('Done.'); % Normal end
end
