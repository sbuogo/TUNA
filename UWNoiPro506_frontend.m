% Front-end for main function 'UWNoiPro' v. 5.06 without GUI
% Setup argument contains all parameters defined by user in GUI
% NOTE: decidecade bands start from 1 Hz, earlier TOB from 10 Hz

Setup.PROCSPL   = false;   % Evaluate SPL (true=yes, false=no)
Setup.PROCPSDB   = true;   % Evaluate PSD (true=yes, false=no)
Setup.PROCPSDD   = true;   % Evaluate TOB (true=yes, false=no)
Setup.BOXPLT    = true;   % TOB plot type: true=boxplot, false=line plot
Setup.STARTTIME = 0;      % Start time of analysis in s
Setup.STOPTIME  = 60;    % Stop time of analysis in s (Inf=end of file)
Setup.SNAPLEN   = 60;     % Snapshot duration in s
Setup.SNAPOVLAP = 0.0;    % Snapshot overlap in percent (e.g. 50.0 = 50%)
Setup.SNAPSKIP  = 0.0;    % Time separation between snapshots in s
Setup.NFFT      = 2^20;   % e.g. 2^10=1024,2^12=4096,2^14=16384,2^16=65536
Setup.FFTOVLAP  = 50;     % FFT overlap in percent: typical 0 or 50
Setup.WTYPE     = 1;      % Window: 0=none, 1=Hamming, 2=Blackman, 3=Hann
Setup.AVGTYPE   = 1;      % Snapshot average: 1=arithmetic mean, 2=median
Setup.HIPASSF   = 0;      % Highpass filter cutoff in Hz (0=no filter)
Setup.PRCTILEL  = 5;      % Lower percentile (e.g. 5 = 5%, 25 = 25%)
Setup.PRCTILEH  = 95;     % Higher percentile (e.g. 75 = 75%, 95 = 95%)
Setup.DDIDXL   = 2;      % First TOB index (see table below)
Setup.DDIDXH   = 21;     % Last TOB index  (see table below)
  %   1=1Hz     11=10Hz    21=100Hz   31=1000Hz   41=10kHz    51=100kHz
  %   2=1.25 Hz 12=12.5Hz  22=125Hz   32=1250Hz   42=12.5kHz  52=125kHz
  %   ...       13=16Hz    23=160Hz   33=1600Hz   43=16kHz    ...
  %   ...       14=20Hz    24=200Hz   34=2000Hz   44=20kHz    ...
  %   ...       15=25Hz    25=250Hz   35=2500Hz   45=25kHz
  %             16=31.5Hz  26=315Hz   36=3150Hz   46=31.5kHz
  %             17=40Hz    27=400Hz   37=4000Hz   47=40kHz
  %             18=50Hz    28=500Hz   38=5000Hz   48=50kHz
  %             19=63Hz    29=630Hz   39=6300Hz   49=63kHz
  %   10=8 Hz   20=80Hz    30=800Hz   40=8000Hz   50=80kHz    60=800kHz
Setup.DRAWALL   = 1;      % Draw all PSD curves (1=yes, 0=only low/avg/hi)
Setup.SAVEJPG   = true;   % Save figure as *.jpg (true=yes,false=no)
Setup.SAVETXT   = true;  % Save results as text file (true=yes,false=no)
Setup.XAXAUTO   = true;   % PSD,TOB X axis autorange (true=yes,false=no)
Setup.YAXAUTO   = false;   % PSD,TOB Y axis autorange (true=yes,false=no)
Setup.XAXMIN    = 5;      % if not autorange, PSD X axis low limit in Hz
Setup.XAXMAX    = 24000;  % if not autorange, PSD X axis high limit in Hz 
Setup.YAXMIN    = 20;     % if not autorange, PSD Y axis low limit in dB
Setup.YAXMAX    = 140;    % if not autorange, PSD Y axis high limit in dB 
Setup.LEGND     = true;  % Show legend on PSD,TOB plot (true=yes,false=no)
Setup.FIGSIZE   = [640 480]; % Size in pixels of all figures
Setup.TITL      = true; % show title in all figures
Setup.VERB      = 2;      % Display info (0=normal,1=detailed,2=debug)
Setup.PLTPAUSE  = 0;    % Pause in s after each plot
Setup.KNDSDRAW  = true;  % true=draw Knudsen curves, false=do not draw
Setup.KNDSSIDX  = [0 1 2];  % sea states of Knudsen curves, 0 to 6
Setup.LOGFILE   = true; % save command output to file
Setup.LOGSETUP  = true; % add setup in log file
Setup.NORESIZ   = true; % keep current figure size and position

Setup.INPATH    = '/home/silvano/ownCloud_lab/UWNoise/test_in'; 
Setup.OUTPATH   = '/home/silvano/ownCloud_lab/UWNoise/test_out';
Setup.CALFPATH='/home/silvano/ownCloud_lab/UWNoise/cal/cal_fix-190db.txt';

[retVal,exStat,figHandles] = UWNoiPro506main(Setup);

disp(exStat)
