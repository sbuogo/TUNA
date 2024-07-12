% evalPSDBBDD() - Subroutine to evaluate Power Spectral Density (PSD)
% either in broadband (BB) or in decidecade bands (DD). 
% Use 'periodogram' with given FFT no. within a single snapshot, as
% many times as needed to fill each snapshot with given overlap.
% Use the given calibration data converted to Scale Factor in Pa.
% Output is normalized to window total power, already accounted for 
% in function 'periodogram()', and to frequency bin width.
% Input arguments: snapshot, FFT window, calib.SF (BB and DD), 
% DD frequencies, parameters. DD arguments are void if not required.
% Output are PSD matrixes in dB re 1 microPa^2/Hz with hi/avg/low values.
% Averaging on snapshots may be either arithm.mean or median.
% Use zero padding for last FFT subinterval, centered and with linear 
% interpolation to remove jump discontinuities.
% DD matrix provides band limits and indexes of f.vector.
% Draw BB or DD plot only if corresponding option is enabled in main.
% Debug mode requires key strokes between plots.
%
% Reference: ISO 18405:2017 3.1.3.13, 3.1.4.2; Lurton p.132
%
% Revision history: see previous 'evalPSD' and 'evalTOB' subroutines.
%
% Current version: 
% 3.0 - 27-mar-2024 - Merged from 'evalPSD26' and 'evalTOB27' subroutines.
%       Changed frequency intervals from one-third octave to decidecade.
%       BB elements are summed within DD band and divided by bandwidth,
%       then normalized to f bin width using srate/nfft factor.
%       Earlier subroutines used band averaging without normalization.
%       Added WAV figure handle and boolean proc.type as input param.
%       Input arguments are void if processing is not required.
%       Two output arguments - BB and DD matrixes - instead of one.
%       Requires main calling function 'UWNoiPro506main' or later.
%
% SBU@INM

function [bbpxxavg,ddpxxavg] = ...
                evalPSDBBDD30(Snap,Win,BBSF,DDSF,DDmat,Param)

% Parameters:
% - Snap.Y   = snapshot array of samples
% - Snap.T   = snapshot time array in s
% - Win      = time window
% - SF       = calibration (scale factor in Pa) for BB frequencies
% - DDSF     = as above, for DD center frequencies
% - DDmat    = decidecade band matrix with center and edge freq.s
% - Param.sr = sample rate
% - Param.nf = FFT n. of points
% - Param.av = averaging type: 1=arithm.mean, 2=median
% - Param.fn = n. of FFTs in each snapshot
% - Param.of = offset in samples between successive FFTs
% - Param.zp = zero pad samples added for last FFT
% - Param.sn = current snapshot number
% - Param.re = pressure reference (normally 1 microPa)
% - Param.vb = verbose mode (0=silent, 1=verbose)
% - Param.pp = pause after each plot in s
% - Param.figwav = figure handle for WAV plot, used in debug mode
% - Param.figpsdb = figure handle for broadband PSD plot
% - Param.figpsdd = figure handle for decidecade PSD plot
% - Param.procpsdb = broadband PSD processing selected (boolean)
% - Param.procpsdd = decidecade PSD processing selected (boolean)

fftno = ceil(Param.fn); % n. of FFTs, including last one with zero pad
PREF = Param.re;

VERB = Param.vb;
PLTPAUSE = Param.pp;
PROCPSDB = Param.procpsdb; % broadband
PROCPSDD = Param.procpsdd; % decidecade

% plot titles and labels
CURTIT = 'Snap %d, FFT n.%d - high/current/low PSD'; % 2 arguments
AVGTIT = 'Snap %d, %d FFTs - high/average/low PSD'; % 2 arguments
BBXLAB = 'Frequency / Hz';
BBYLAB = 'PSD / dB re 1 \muPa^2/Hz';
DDXLAB = BBXLAB;
DDYLAB = BBYLAB;

figWAV = Param.figwav; % WAV figure handle
if PROCPSDB
    figPSDB = Param.figpsdb; % broadband PSD figure handle
end
if PROCPSDD
    figPSDD = Param.figpsdd; % decidecade PSD figure handle
end

DOTN = 10; % display '.' character every N cycles to show progress

fbinnorm = (Param.sr/Param.nf); % normalize for f.bin width (linear)
zeropadnorm = (1 - Param.zp/Param.nf); % normalize for zero pad

if VERB > 0    % add details if verbose mode on
    fprintf('\nParam: sr=%g, nf=%d, fn=%g, of=%g, zp=%d', ...
        Param.sr, ...
        Param.nf, ...
        Param.fn, ...
        Param.of, ...
        Param.zp)
    fprintf('\nNorm: fbin=%g, zp=%g',fbinnorm,zeropadnorm); 
end

bbpxxmat = zeros(Param.nf/2 + 1, fftno); % allocate BB buffer matrix

if PROCPSDD
    ddf = DDmat(:,1); % vector of TOB center freq.s, used in plots
    ddpxxmat = zeros(size(DDmat,1), fftno); % allocate TOB buffer
end

%%%%%%%%%%%%%%%%% loop over n. of FFTs within snapshot %%%%%%%%%%%%%%%%
for i = 1:fftno % last portion might be fractional
    
    startSmp = 1 + (i-1) * Param.of;
    stopSmp  = startSmp + Param.nf - 1;
    
    if stopSmp <= length(Snap.Y) % normal FFT vector within snapshot
        Y = Snap.Y(startSmp:stopSmp); 
        zpString = '';
    else % last FFT: use zero padding
        % done adding half zeros before and half after, to smooth
        % possible start/endpoint discontinuities by windowing
        % note: zeropadding should be with even n. of points
        zplen = Param.zp / 2;
        zPad = zeros(zplen, 1);
        % linearly interpolate in zeropad subarrays
        rampUp = interp1([1,zplen],[0,Snap.Y(startSmp)],1:zplen);
        rampDn = interp1([1,zplen],[Snap.Y(length(Snap.Y)),0],1:zplen);
        Ynointerp = [zPad; Snap.Y(startSmp:length(Snap.Y)); zPad];
        Y = [rampUp'; Snap.Y(startSmp:length(Snap.Y)); rampDn'];
        zpString = sprintf('zp=%d',Param.zp);
    end
        
    % DEBUG: draw Y (with and without zeropad interpolation), hit any key
    if VERB == 2
        set(0,'CurrentFigure',figWAV); % draw on waveform figure
        if length(zpString) > 0
            plot(1:Param.nf,Ynointerp)
            pause
        end
        plot(1:Param.nf,Y)
        pause
    end
    % END DEBUG
    
    if VERB > 0 % display progress info
        fprintf('\n%d: [%d %d] %s',i,startSmp,stopSmp,zpString); 
    else
        if mod(i,DOTN) == 0
            fprintf('.'); % silent progress, only print dot every DOTN
        end
    end

    %%%%%%%%%%%%% compute raw PSD through periodogram on windowed data;
    % output is normalized to power of window specified as parameter

    [pxx,bbf] = periodogram(Y, Win, Param.nf, Param.sr);
    
    if i == fftno && Param.zp > 0 % last FFT with zero padding
        pxx = pxx / zeropadnorm; % normalize to zeropad length
    end

    bbpxxSF = (pxx .* (BBSF'.^2)); % apply BB calibration scale factor

    bbpxxmat(:,i) = bbpxxSF; % update column vector in buffer matrix

    %%%%%%%% BB min-max bounds, to be repeated over DDecs if enabled
    bbpxxdb = 10*log10(bbpxxSF) - 20*log10(PREF); % convert to dB
    if VERB > 0
        fprintf(' max = %g dB',max(bbpxxdb)); % disp. max if verbose
    end

    if i == 1   % initialize max & min
        bbpxxdbhi = bbpxxdb;
        bbpxxdblo = bbpxxdb;
    else           % update max & min
        for n = 1 : length(bbpxxdb)  
            if bbpxxdb(n) > bbpxxdbhi(n)
                bbpxxdbhi(n) = bbpxxdb(n);
            end
            if bbpxxdb(n) < bbpxxdblo(n)
                bbpxxdblo(n) = bbpxxdb(n);
            end
        end
    end % BB min-max bounds

    if PROCPSDD %%%%%%%%% decidecade min-max bounds, repeat from above.
        % For each DD band, add BB elements & divide by DD bandwidth.
        for ddidx = 1:size(DDmat,1) % DD index
            fvil = int32(DDmat(ddidx,6)); % lower DD limit f.vect.idx
            fvih = int32(DDmat(ddidx,7)); % upper DD limit f.vect.idx
            ddlow = DDmat(ddidx,2);       % nominal DD lower limit in Hz
            ddhi  = DDmat(ddidx,3);       % nominal DD upper limit in Hz
            
            % f bin normalization factor s.rate/n.fft: if factor=1
            % then each element is exactly 1 Hz wide, otherwise need to
            % account for bin width to compute total power.
            ddpxx(ddidx,i) = sum(pxx(fvil:fvih))*fbinnorm/(ddhi-ddlow);
            
            % OLD     ddpxx(ddidx,i) = mean(pxx(fvil:fvih));
            % in the old version above, taking the mean was equivalent to
            % summing and dividing by (b.width in Hz / bin width in Hz)
            % where the quantity in () is the n. of DD band elements
        end

        ddpxxSF(:,i) = ddpxx(:,i) .* (DDSF.^2); % decidecade calib.SF
        ddpxxdb = 10*log10(ddpxxSF(:,i)) - 20*log10(PREF); % to dB
        
        if i == 1   % initialize max & min on first FFT cycle
            ddpxxdbhi = ddpxxdb;
            ddpxxdblo = ddpxxdb;
        else
            for n = 1 : size(ddpxxmat,1)   % update max & min
                if ddpxxdb(n) > ddpxxdbhi(n)
                    ddpxxdbhi(n) = ddpxxdb(n);
                end
                if ddpxxdb(n) < ddpxxdblo(n)
                    ddpxxdblo(n) = ddpxxdb(n);
                end
            end
        end
    end % decidecade min-max bounds

    %%%%%% plot hi-current-low only on enabled figures,
         % although BB data are always present
    if PROCPSDB
        set(0,'CurrentFigure',figPSDB); % draw on BB PSD figure
        semilogx(bbf,bbpxxdbhi,'r',bbf,bbpxxdblo,'g',bbf,bbpxxdb,'b' )
        grid on
        title(sprintf(CURTIT,Param.sn,i));
        xlabel(BBXLAB)
        ylabel(BBYLAB)
    end
    if PROCPSDD
        set(0,'CurrentFigure',figPSDD); % draw on DD PSD figure
        semilogx(ddf,ddpxxdbhi,':or', ...
                 ddf,ddpxxdblo,':og', ...
                 ddf,ddpxxdb,'ob' )
        grid on
        title(sprintf(CURTIT,Param.sn,i));
        xlabel(DDXLAB)
        ylabel(DDYLAB)
    end
    
    if VERB == 2
        pause %%%%%% DEBUG: hit any key to continue
    else
        pause(PLTPAUSE/10) % normal run should be fast, short pause
    end

end %%%%%%%%%%%%%% FFT loop %%%%%%%%%%%%%%%%%

%%%%%% plot hi-low-averages for BB, TOB, if enabled

if PROCPSDB % broadband
    switch Param.av % average on linear data
        case 1 % arithmetic mean
            bbpxxavg =   mean(bbpxxSF,2);
        case 2 % median
            bbpxxavg = median(bbpxxSF,2);
    end
    bbpxxdbavg = 10*log10(bbpxxavg) - 20*log10(PREF); 
    set(0,'CurrentFigure',figPSDB);
    semilogx(bbf,bbpxxdbhi,'r', bbf,bbpxxdblo,'g',bbf,bbpxxdbavg,'b')
    grid on
    title(sprintf(AVGTIT,Param.sn,i));
else
    % BB processing done anyway, should not need to return void
end

if PROCPSDD
    switch Param.av % average on linear data
        case 1 % arithmetic mean
            ddpxxavg =   mean(ddpxxSF,2);
        case 2 % median
            ddpxxavg = median(ddpxxSF,2);
    end
    ddpxxdbavg = 10*log10(ddpxxavg) - 20*log10(PREF); 
    set(0,'CurrentFigure',figPSDD);
    semilogx(ddf,ddpxxdbhi,':or',ddf,ddpxxdblo,':og') % plot hi-low
    hold on
    semilogx(ddf,ddpxxdbavg,'ob','MarkerFaceColor','b') % plot avg
    hold off
    grid on
    title(sprintf(AVGTIT,Param.sn,i));
else
    ddpxxavg = []; % return void if TOB option not selected in main
end

if VERB == 2
    pause % DEBUG: hit any key to continue
else
    pause(PLTPAUSE) % normal run, pause as user specified
end
return
