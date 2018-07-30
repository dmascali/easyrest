function varargout = cpsd_R2014b(x,y,varargin)
%CPSD   Cross Power Spectral Density (CPSD) estimate via Welch's method.
%   Pxy = CPSD(X,Y) returns the Cross Power Spectral Density estimate, Pxy,
%   of two discrete-time signals, X and Y, using Welch's averaged,
%   modified periodogram method. By default, X and Y are divided into
%   eight sections with 50% overlap, each section is windowed with a
%   Hamming window and eight modified periodograms are computed and
%   averaged. See "help pwelch" and "help cpsd" for complete details.
%
%   X and Y may be either vectors or two-dimensional matrices. If both are
%   matrices, they must have the same size, and CPSD operates columnwise:
%   Pxy(:,n) = CPSD(X(:,n),Y(:,n)). If one is a matrix and the other is a
%   vector, the vector is converted to a column vector and internally
%   expanded so both inputs have the same number of columns.
%
%   Pxy is the distribution of power per unit frequency. For real signals,
%   CPSD returns the one-sided Cross PSD by default; for complex signals,
%   it returns the two-sided Cross PSD. Note that a one-sided Cross PSD
%   contains the total power of the input signal.
%
%   Pxy = CPSD(X,Y,WINDOW), when WINDOW is a vector, divides each column of
%   X and Y into overlapping sections of length equal to the length of
%   WINDOW, and then windows each section with the vector specified in
%   WINDOW. If WINDOW is an integer, then each column of X and Y are
%   divided into sections of length WINDOW, and each section is windowed
%   with a Hamming of that length. If WINDOW is omitted or specified as
%   empty, a Hamming window is used to obtain eight sections of X and Y.
%
%   Pxy = CPSD(X,Y,WINDOW,NOVERLAP) uses NOVERLAP samples of overlap from
%   section to section. NOVERLAP must be an integer smaller than the
%   length of WINDOW if WINDOW is a vector, or smaller than WINDOW if
%   WINDOW is an integer. If NOVERLAP is omitted or specified as empty, it
%   is set to obtain a 50% overlap.
%
%   [Pxy,W] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT) specifies the number of FFT
%   points used to calculate the Cross PSD estimate. For real signals, Pxy
%   has length (NFFT/2+1) if NFFT is even, and (NFFT+1)/2 if NFFT is odd.
%   For complex signals, Pxy always has length NFFT. If NFFT is specified
%   as empty, the  default NFFT -the maximum of 256 or the next power of
%   two greater than the length of each section of X (and Y)- is used.
%
%   If NFFT is greater than the length of each section, the data is
%   zero-padded. If NFFT is less than the section length, the segment is
%   "wrapped" (using DATAWRAP) to make the length equal to NFFT. This
%   produces the correct FFT when NFFT is smaller than the section length.
%
%   W is the vector of normalized frequencies at which the PSD is
%   estimated. W has units of radians/sample. For real signals, W spans
%   the interval [0,pi] when NFFT is even and [0,pi) when NFFT is odd. For
%   complex signals, W always spans the interval [0,2*pi).
%
%   [Pxy,W] = CPSD(X,Y,WINDOW,NOVERLAP,W) uses the Goertzel algorithm to
%   compute the two-sided Cross CPSD at the normalized angular frequencies
%   contained in the vector W. W must have at least two elements. The
%   frequencies specified in W are rounded to the nearest DFT bin
%   commensurate with the signal's resolution.
%
%   [Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs) returns the Cross PSD as a
%   function of physical frequency. Fs is the sampling frequency specified
%   in hertz. If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which Pxy is estimated.
%   For real signals, F spans the interval [0,Fs/2] when NFFT is even and
%   [0,Fs/2) when NFFT is odd. For complex signals, F always spans the
%   interval [0,Fs).
%
%   [Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,F,Fs) uses the Goertzel algorithm to
%   compute the Cross PSD estimate at the physical frequencies contained in
%   the vector F. F must be expressed in hertz and have at least two
%   elements. The frequencies specified in F are rounded to the nearest DFT
%   bin commensurate with the signal's resolution.
%
%   [...] = CPSD(...,FREQRANGE) returns the Cross PSD computed over the
%   specified range of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided Cross PSD of real input signals X
%         and Y. If NFFT is even, Pxy has length NFFT/2+1 and is computed
%         over the interval [0,pi]. If NFFT is odd, Pxy has length
%         (NFFT+1)/2 and is computed over the interval [0,pi). When Fs is
%         optionally specified, the intervals become [0,Fs/2) and [0,Fs/2]
%         for even and odd NFFT, respectively.
%
%      'twosided' - returns the two-sided Cross PSD for either real or
%         complex input X and Y. Pxy has length NFFT and is computed over
%         the interval [0,2*pi). When Fs is specified, the interval becomes
%         [0,Fs).
%
%      'centered' - returns the centered two-sided Cross PSD for either
%         real or complex X and Y. Pxy has length NFFT and is computed
%         over the interval (-pi, pi] for even NFFT and (-pi, pi) for odd
%         NFFT. When Fs is specified, the intervals become (-Fs/2, Fs/2]
%         and (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP. The default value of FREQRANGE is 'onesided' when X
%      and Y are both real and 'twosided' when either X or Y is complex.
%
%   CPSD(...) with no output arguments plots the Cross PSD (in decibels per
%   unit frequency) in the current figure window.
%
%   EXAMPLE:
%      Fs = 1000;   t = 0:1/Fs:.296;
%      x = cos(2*pi*t*200)+randn(size(t));  % A cosine of 200Hz plus noise
%      y = cos(2*pi*t*100)+randn(size(t));  % A cosine of 100Hz plus noise
%      cpsd(x,y,[],[],[],Fs,'twosided');    % Uses default window, overlap & NFFT. 
% 
%   See also PWELCH, PERIODOGRAM, PCOV, PMCOV, PBURG, PYULEAR, PEIG, PMTM,
%   PMUSIC.

%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and 
%         Modeling, John Wiley & Sons, 1996.

narginchk(1,7);
nargoutchk(0,3);

esttype = 'cpsd';
% Possible outputs are:
%       Plot
%       Pxx
%       Pxx, freq
[varargout{1:nargout}] = welch({x,y},esttype,varargin{:});

% [EOF]


function varargout = computepsd(Sxx,w,range,nfft,Fs,esttype)
%COMPUTEPSD  Compute the one-sided or two-sided PSD or Mean-Square.
%   [Pxx,W,UNITS] = COMPUTEPSD(Sxx,W,RANGE,NFFT,Fs,ESTTYPE) where the
%   inputs and outputs are:
%
%   Inputs:
%    Sxx   - Whole power spectrum [Power]; it can be a vector or a matrix.
%            For matrices the operation is applied to each column.
%    W     - Frequency vector in rad/sample or in Hz.
%    RANGE - Determines if a 'onesided' or a 'twosided' Pxx and Sxx are
%            returned.
%    NFFT  - Number of frequency points.
%    Fs    - Sampling Frequency.
%    ESTTYPE - A string indicating the estimate type: 'psd', or 'ms' value.
%
%   Outputs:
%    Pxx   - One-sided or two-sided PSD or MEAN-SQUARE (not scaled by Fs)
%            depending on the input arguments RANGE and TYPE.
%    W     - Frequency vector 0 to 2*Nyquist or 0 to Nyquist depending on
%            range, units will be either rad/sample (if Fs is empty) or Hz
%            (otherwise).
%    UNITS - Either 'rad/sample' or 'Hz' 

%   Author(s): R. Losada
%   Copyright 1988-2012 The MathWorks, Inc.

if nargin < 6,
    esttype = 'psd';
end

w = w(:); % Make sure we always returns a column vector for frequency

% Generate the one-sided spectrum [Power] if so wanted
if strcmp(range,'onesided'),
   if rem(nfft,2),
      select = 1:(nfft+1)/2;  % ODD
      Sxx_unscaled = Sxx(select,:); % Take only [0,pi] or [0,pi)
      Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled
   else
      select = 1:nfft/2+1;    % EVEN
      Sxx_unscaled = Sxx(select,:); % Take only [0,pi] or [0,pi)
      Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end-1,:); Sxx_unscaled(end,:)]; % Don't double unique Nyquist point
   end
   w = w(select);
end

% Compute the PSD [Power/freq]
if ~isempty(Fs),
   Pxx = Sxx./Fs; % Scale by the sampling frequency to obtain the psd
   units = 'Hz';
else
   Pxx = Sxx./(2.*pi); % Scale the power spectrum by 2*pi to obtain the psd
   units = 'rad/sample';    
end

if any(strcmpi(esttype,{'ms','power'}))
    varargout = {Sxx,w,units};  % Mean-square
else  
    varargout = {Pxx,w,units};  % PSD 
end

% [EOF] computepsd.m

function w = psdfreqvec(varargin)
%This undocumented function may be removed in a future release.
%
%PSDFREQVEC Frequency vector
%   PSDFREQVEC('Npts',NPTS) returns a frequency vector in radians based on
%   the number of points specified in NPTS. The vector returned assumes 2pi
%   periodicity.
%
%   PSDFREQVEC('Fs',FS) specifies the sampling frequency FS in hertz. By
%   default Fs is set to empty indicating normalized frequency.
%   
%   PSDFREQVEC('CenterDC',CENTERDC) specifies a boolean value in CENTERDC
%   which indicates if zero hertz should be in the center of the frequency
%   vector. CENTERDC can be one of these two values [ {false} | true ]. 
%
%   PSDFREQVEC('Range',RANGE) specifies the range of frequency in RANGE.
%   RANGE can be one of the two strings [ {whole} | half ]. Assuming
%   CenterDC=false then:
%       'whole' = [0, 2pi)
%       'half'  = [0, pi] for even NPTS or [0, pi) for odd NPTS
%
%   When CenterDC=true then:
%       'whole' = (-pi, pi] for even NPTs or (-pi, pi) for odd NPTs
%       'half'  = [-pi/2, pi/2] for even* NPTS or (-pi/2, pi/2) for odd NPTS
%
%       *When NPTS is not divisible by 4, then the range is (-pi/2, pi/2).
%
%   When Range='half' the frequency vector has length (NPTS/2+1) if NPTS is
%   even**, and (NPTS+1)/2 if NPTS is odd***.
%
%       **If CenterDc=true and the number of points specified is even is
%       not divisible by 4, then the number of points returned is NPTS/2.
%       This is to avoid frequency points outside the range [-pi/2 pi/2]. 
%
%       ***If CenterDC=true and the number of points NPTS specified is odd
%       and (NPTS+1)/2 is even then the length of the frequency vector is
%       (NPTS-1)/2.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

% error(nargchk(1,3,nargin,'struct'));

% Define default parameters; use same parameter names as defined in help.
Npts = 1024;
Fs = [];
Range = 'whole';
CenterDC = false;
local_pvparse(varargin{:});

% Compute the frequency grid.
w = frequencygrid(Range,Npts,Fs,CenterDC);

%--------------------------------------------------------------------------
function w = frequencygrid(Range,Npts,Fs,CenterDC)
% Compute the frequency grid.

% Compute the whole frequency range, e.g., [0,2pi) to avoid round off errors.
if isempty(Fs),
    Fs = 2*pi;
end
freq_res = Fs/Npts;
w = freq_res*(0:Npts-1);

% There can still be some minor round off errors in the frequency grid.  
% Fix the known points, i.e., those near pi and 2pi.
Nyq = Fs/2;
half_res = freq_res/2; % half the resolution

% Determine if Npts is odd and calculate half and quarter of Npts.
[isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(Npts);

if isNPTSodd,
    % Adjust points on either side of Nyquist.
    w(halfNPTS)   = Nyq - half_res;
    w(halfNPTS+1) = Nyq + half_res;
else
    % Make sure we hit Nyquist exactly, i.e., pi or Fs/2 
    w(halfNPTS) = Nyq;
end
w(Npts) = Fs-freq_res;

% Get the right grid based on range, centerdc, etc.
w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS);

%--------------------------------------------------------------------------
function [isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(NPTS)
% Determine if we're dealing with even or odd lengths of NPTS, 1/2 NPTS,
% and 1/4 NPTS.

% Determine if Npts is odd.
isNPTSodd = false;
if rem(NPTS,2),
    isNPTSodd = true;
end

% Determine half the number of points.
if isNPTSodd,   halfNPTS = (NPTS+1)/2;  % ODD
else            halfNPTS = (NPTS/2)+1;  % EVEN
end

% Determine if half Npts is odd.
ishalfNPTSodd = false;     
if rem(halfNPTS,2),        
    ishalfNPTSodd = true;  
end

% Determine a quarter of the number of points.
if ishalfNPTSodd,  quarterNPTS = (halfNPTS+1)/2;  % ODD
else               quarterNPTS = (halfNPTS/2)+1;  % EVEN
end

%--------------------------------------------------------------------------
function w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS)
% Calculate the correct grid based on user specified values for range,
% centerdc, etc.

switch lower(Range)
    case 'whole',
        % Calculated by default.% [0, 2pi)

        if CenterDC,          % (-pi, pi] even or (-pi, pi) odd
            if isNPTSodd,  negEndPt = halfNPTS;
            else           negEndPt = halfNPTS-1;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:halfNPTS)];
        end
        
    case 'half'            
        w = w(1:halfNPTS);      % [0, pi] even or [0, pi) odd
        
        % For even number of points that are not divisible by 4 you get
        % less one point to avoid going outside the [-pi/2 pi/2] range.
        if CenterDC,            % [-pi/2, pi/2] even (-pi/2, pi/2) odd 
            if ishalfNPTSodd,
                negEndPt = quarterNPTS;
            else
                quarterNPTS = quarterNPTS-1; % Avoid going over pi/2
                negEndPt = quarterNPTS;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:quarterNPTS)];
            if ~rem(Npts,4),
                % Make sure we hit pi/2 exactly when Npts is divisible
                % by 4! In this case it's due to roundoff.
                w(end) = Nyq/2;
            end
        end
    otherwise
        error(message('signal:psdfreqvec:InternalError'));
end
w = w(:);  % Return a column vector.

%--------------------------------------------------------------------------
function local_pvparse(varargin)

varnames = evalin('caller','whos');
vnames = {varnames.name};
for m = 1:2:nargin
    indx = find(strncmpi(vnames, varargin{m}, length(varargin{m})));
    switch length(indx)
        case 0
            error(message('signal:psdfreqvec:unknownInput', varargin{ m }));
        case 1
            assignin('caller', vnames{indx}, varargin{m+1});
        otherwise
            error(message('signal:psdfreqvec:ambiguousInput', varargin{ m }));
    end
end

% [EOF]

function [Xx,f] = computeDFT(xin,nfft,varargin)
%COMPUTEDFT Computes DFT using FFT or Goertzel
%   This function is used to calculate the DFT of a signal using the FFT 
%   or the Goertzel algorithm. 
%
%   [XX,F] = COMPUTEDFT(XIN,NFFT) where NFFT is a scalar and computes the 
%   DFT XX using FFT. F is the frequency points at which the XX is 
%   computed and is of length NFFT.
%
%   [XX,F] = COMPUTEDFT(XIN,F) where F is a vector with atleast two 
%   elements computes the DFT XX using the Goertzel algorithm. 
%
%   [XX,F] = COMPUTEDFT(...,Fs) returns the frequency vector F (in hz)
%   where Fs is the sampling frequency
%
%   Inputs:
%   XIN is the input signal
%   NFFT if a scalar corresponds to the number of FFT points used to 
%   calculate the DFT using FFT.
%   NFFT if a vector corresponds to the frequency points at which the DFT
%   is calculated using goertzel.
%   FS is the sampling frequency 

% Copyright 2006-2014 The MathWorks, Inc.

% [1] Oppenheim, A.V., and R.W. Schafer, Discrete-Time Signal Processing,
% Prentice-Hall, Englewood Cliffs, NJ, 1989, pp. 713-718.
% [2] Mitra, S. K., Digital Signal Processing. A Computer-Based Approach.
% 2nd Ed. McGraw-Hill, N.Y., 2001.

narginchk(2,3);
if nargin > 2,
    Fs = varargin{1};
else
    Fs = 2*pi;
end

nx = size(xin,1);

if length(nfft) > 1, 
    isfreqVector = true; 
else 
    isfreqVector = false;
end

if ~isfreqVector,
    [Xx,f] = computeDFTviaFFT(xin,nx,nfft,Fs);
else
    [Xx,f] = computeDFTviaGoertzel(xin,nfft,Fs);
end
    
%-------------------------------------------------------------------------
function [Xx,f] = computeDFTviaFFT(xin,nx,nfft,Fs)
% Use FFT to compute raw STFT and return the F vector.

% Handle the case where NFFT is less than the segment length, i.e., "wrap"
% the data as appropriate.
xin_ncol = size(xin,2);
xw = zeros(nfft,xin_ncol);
if nx > nfft,
    for j = 1:xin_ncol, 
        xw(:,j) = datawrap(xin(:,j),nfft);
    end
else
    xw = xin;
end

Xx = fft(xw,nfft);
f = psdfreqvec('npts',nfft,'Fs',Fs);
%--------------------------------------------------------------------------
function [Xx,f] = computeDFTviaGoertzel(xin,freqvec,Fs)
% Use Goertzel to compute raw DFT and return the F vector.

f = freqvec(:);
f = mod(f,Fs); % 0 <= f < = Fs
xm = size(xin,1); % NFFT

% Indices used by the Goertzel function (see equation 11.1 pg. 755 of [2])
fscaled = f/Fs*xm+1;
k = round(fscaled);

% shift for each frequency from default xm length grid
deltak = fscaled-k;

tempk = k;
% If k > xm, fold over to the 1st bin
k(tempk > xm) = 1;

n = (0:xm-1)';
Xx = zeros(size(k,1),size(xin,2));
for kindex = 1:length(k)
    % We need to evaluate the DFT at the requested frequency instead of a
    % neighboring frequency that lies on the grid obtained with xm number
    % of points in the 0 to fs range. We do that by giving a complex phase
    % to xin equal to the offset from the frequency to its nearest neighbor
    % on the grid. This phase translates into a shift in the DFT by the
    % same amount. The Xx(k) then is the DFT at (k+deltak).
    
    % apply kernal to xin so as to evaluate DFT at k+deltak)
    kernel = exp(-1i*2*pi*deltak(kindex)/xm*n);
    xin_phaseshifted = xin.*repmat(kernel,1,size(xin,2));
    
    Xx(kindex,:) = goertzel(xin_phaseshifted,k(kindex));
end

% DFT computed at exactly the frequencies it was requested for
f = freqvec(:);

%-----------------------------

function [P,f] = computeperiodogram(x,win,nfft,esttype,varargin)
%COMPUTEPERIODOGRAM   Periodogram spectral estimation.
%   This function is used to calculate the Power Spectrum Sxx, and the
%   Cross Power Spectrum Sxy.
%
%   Sxx = COMPUTEPERIODOGRAM(X,WIN,NFFT) where x is a vector returns the
%   Power Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Sxy = COMPUTEPERIODOGRAM({X,Y},WIN,NFFT) returns the Cross Power
%   Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Inputs:
%    X           - Signal vector or a cell array of two elements containing
%                  two signal vectors.
%    WIN         - Window
%    NFFT        - Number of frequency points (FFT) or vector of
%    frequencies at which periodogram is desired (Goertzel)
%    WINCOMPFLAG - A string indicating the type of window compensation to
%                  be done. The choices are: 
%                  'ms'    - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'power' - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'psd'   - compensate for Power Spectral Density (PSD);
%                            maintain correct area under the PSD curve.
%
%   Output:
%    Sxx         - Power spectrum [Power] over the whole Nyquist interval. 
%      or
%    Sxy         - Cross power spectrum [Power] over the whole Nyquist
%                  interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(3,5);
if nargin < 4,
    esttype = 'psd'; % Default, compensate for window's power.
end

if nargin < 5 || isempty(varargin{1}),
    Fs = 2*pi;
else
    Fs = varargin{1};     
end

% Validate inputs and convert row vectors to column vectors.
[x,~,y,is2sig,win] = validateinputs(x,win,nfft);

% Window the data
xw = bsxfun(@times,x,win);
if is2sig, yw = bsxfun(@times,y,win); end 

% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if any(strcmpi(esttype,{'ms','power'}))
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win)^2;
else
    U = win'*win;  % compensates for the power of the window.
end

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels

[Xx,f] = computeDFT(xw,nfft,Fs);
if is2sig, [Yy,f] = computeDFT(yw,nfft,Fs); end

P = Xx.*conj(Xx)/U;      % Auto spectrum.
if is2sig,
    P = bsxfun(@times,Xx,conj(Yy))/U;  % Cross spectrum.
end

%--------------------------------------------------------------------------
function [x,Lx,y,is2sig,win] = validateinputs(x,win,~)
% Validate the inputs to computexperiodogram and convert row vectors to
% column vectors for backwards compatiblility with R2014a and prior
% releases

% Set defaults and convert to row vectors to columns.
y     = [];
is2sig= false;
win   = win(:);
Lw    = length(win);

% Determine if one or two signal vectors was specified.
if iscell(x),
    if length(x) > 1,
        y = x{2};
        if isvector(y)
            y = y(:);
        end
        is2sig = true;
    end
    x = x{1};
end

if isvector(x)
    x = x(:);
end

Lx = size(x,1);

if is2sig,
    Ly  = size(y,1);
    if Lx ~= Ly,
        error(message('signal:computeperiodogram:invalidInputSignalLength'))
    end
    if size(x,2)~=1 && size(y,2)~=1 && size(x,2) ~= size(y,2)
        error(message('signal:computeperiodogram:MismatchedNumberOfChannels'))
    end
end

if Lx ~= Lw,
    error(message('signal:computeperiodogram:invalidWindow', 'WINDOW'))
end

if (numel(x)<2 || numel(size(x))>2)
    error(message('signal:computeperiodogram:NDMatrixUnsupported'))
end

function flag = sigcheckfloattype(x, dataType, fcnName, varName, datacheckflag)
%SIGCHECKFLOATTYPE Check if input x is floating point or numeric and of the
%expected data type
%
% Inputs:
% x             - input data
% dataType      - data type we want to check ('single','double','int8',...)
%                 if set to empty ('') then we do not check for a specific
%                 data type. We only check if data is floating point or
%                 numeric depending on the datacheckflag input.
% fcnName       - function name
% varName       - variable name
% datacheckflag - can be 'allowfloat' or 'allownumeric'. Default is
%                 'allowfloat'. When set to 'allowfloat' the function
%                 checks if data is floating point and then checks if data
%                 is of the specified dataType type. When set to
%                 'allownumeric' the function checks if data is numeric and
%                 then checks if data is of the specified dataType type.
%
% Outputs:
% flag          - true if data is of type dataType


%   Copyright 2013 The MathWorks, Inc.

if nargin < 3
  fcnName = '';
  varName = '';
  datacheckflag = 'allowfloat';
end
if  nargin < 4
  varName = '';
  datacheckflag = 'allowfloat';
end
if nargin < 5
  datacheckflag = 'allowfloat';
end

if strcmpi(datacheckflag,'allowfloat')
  typeCheck = isfloat(x);
  expType = 'double/single';
elseif strcmpi(datacheckflag,'allownumeric')
  typeCheck = isnumeric(x);
  expType = 'numeric';
else
  error(message('signal:sigcheckfloattype:InvalidDataCheckFlag'));    
end

if ~typeCheck
  if ~isempty(fcnName)
    if ~isempty(varName)
      error(message('signal:sigcheckfloattype:InvalidInput',...
        varName, fcnName, expType, class(x)));    
    else
      error(message('signal:sigcheckfloattype:InvalidInput1',...
        fcnName, expType, class(x)));    
    end
  else
    if ~isempty(varName)
      error(message('signal:sigcheckfloattype:InvalidInput2',...
        varName, expType, class(x)));
    else
      error(message('signal:sigcheckfloattype:InvalidInput3',...
        expType, class(x)));
    end
  end    
end

flag = isa(x,dataType);


function y = sigcasttofloat(x, castType, fcnName, varName,datacheckflag)
%SIGCASTTOFLOAT Check if input x is floating point or numeric and then
%casts input to castType
% Inputs:
% x             - input data
% castType      - data type we want to cast to ('single','double')                 
% fcnName       - function name
% varName       - variable name
% datacheckflag - can be 'allowfloat' or 'allownumeric'. Default is
%                 'allowfloat'. When set to 'allowfloat' the function
%                 checks if data is floating point and then casts data to
%                 the specified castType type. When set to 'allownumeric'
%                 the function checks if data is numeric and then casts
%                 data to the specified castType type.
%
% Outputs:
% y            - cast output data

%   Copyright 2013 The MathWorks, Inc.

if nargin < 3
  fcnName = '';
  varName = '';
  datacheckflag = 'allowfloat';
end
if  nargin < 4
  varName = '';
  datacheckflag = 'allowfloat';
end
if nargin < 5
  datacheckflag = 'allowfloat';
end

sigcheckfloattype(x,'', fcnName, varName, datacheckflag);

y = cast(x,castType);

function [options,msg,msgobj] = psdoptions(isreal_x,options,varargin)
%PSDOPTIONS   Parse the optional inputs to most psd related functions.
%
%
%  Inputs:
%   isreal_x             - flag indicating if the signal is real or complex
%   options              - the same structure that will be returned with the
%                          fields set to the default values
%   varargin             - optional input arguments to the calling function
%
%  Outputs:
%  PSDOPTIONS returns a structure, OPTIONS, with following fields:
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided'
%   options.centerdc     - true if 'centered' specified
%   options.conflevel    - confidence level between 0 and 1 ('omitted' when unspecified)
%   options.ConfInt      - this field only exists when called by PMTM 
%   options.MTMethod     - this field only exists when called by PMTM
%   options.NW           - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.Noverlap     - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.CorrFlag     - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.EVFlag       - this field only exists when called by PMUSIC (see PMUSIC for explanation)



%  Authors: D. Orofino and R. Losada
%  Copyright 1988-2012 The MathWorks, Inc.
   
[s,msg,msgobj] = psd_parse(varargin{:});
if ~isempty(msg),
   return
end

omit = 'omitted';

% Replace defaults when appropriate

% If empty or omitted, nfft = default nfft (default is calling function dependent)
if ~any([isempty(s.NFFT),strcmp(s.NFFT,omit)]),    
   options.nfft = s.NFFT;   
end

if ~strcmp(s.Fs,omit), % If omitted, Fs = [], work in rad/sample
   options.Fs = s.Fs;
   if isempty(options.Fs),
      options.Fs = 1;  % If Fs specified as [], use 1 Hz as default
   end
end

% Only PMTM has the field options.ConfInt; 
if isfield(options,'ConfInt'),
   if ~strcmp(s.ConfInt,omit) && ~isempty(s.ConfInt),
      % Override the default option ONLY if user specified a new value
      options.ConfInt = s.ConfInt;
   end
else
   % If a third scalar was specified, error out unless PMUSIC is the calling function
   % (hence options.NW exists)
   if ~strcmp(s.ConfInt,omit) && ~isfield(options,'nw'),  
      msgobj = message('signal:psdoptions:TooManyNumericOptions');
      msg = getString(msgobj);
      return
   end
end

% Only PMUSIC has the field options.nw; 
if isfield(options,'nw'),
   if ~strcmp(s.nw,omit) && ~isempty(s.nw),
      % Override the default option ONLY if user specified a new value
      options.nw = s.nw;
      if ~any(size(options.nw)==1),
         msgobj = message('signal:psdoptions:MustBeScalarOrVector','NW');
         msg = getString(msgobj);
         return
      elseif length(options.nw) > 1,
         options.window = options.nw;
         options.nw = length(options.nw);
      end
   end
else
   % If a third scalar was specified, error out unless PMTM is the calling function
   % (hence options.ConfInt exists)
   if ~strcmp(s.nw,omit) && ~isfield(options,'ConfInt'),  
      msgobj = message('signal:psdoptions:TooManyNumericOptions');
      msg = getString(msgobj);
      return
   end
end

% Only PMUSIC has the field options.noverlap; 
if isfield(options,'noverlap'),
   if ~strcmp(s.noverlap,omit) && ~isempty(s.noverlap),
      % Override the default option ONLY if user specified a new value
      options.noverlap = s.noverlap;
   else
      % Use default
      options.noverlap = options.nw -1;
   end
end
 
options.centerdc = ~strcmp(s.DCFlag,'omitted');
options.conflevel = s.ConfLevel;

if ~strcmp(s.Range,omit),
   options.range = s.Range;
end

if ~isreal_x & strcmpi(options.range,'onesided'), %#ok
   msgobj = message('signal:psdoptions:ComplexInputDoesNotHaveOnesidedPSD');
   msg = getString(msgobj);
   return
end

% Only PMTM has the field options.MTMethod
if ~isfield(options,'MTMethod'),
   if ~strcmp(s.MTMethod,omit),
      % A string particular to pmtm is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.MTMethod,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string
   options.MTMethod = s.MTMethod; 
end

% Only PMUSIC has the field options.CorrFlag
if ~isfield(options,'CorrFlag'),
   if ~strcmp(s.CorrFlag,omit),
      % A string particular to pmusic is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.CorrFlag,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string, we set the flag to one
   options.CorrFlag = 1; 
end

% Only PMUSIC has the field options.EVFlag
if ~isfield(options,'EVFlag'),
   if ~strcmp(s.EVFlag,omit),
      % A string particular to pmusic is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.EVFlag,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string, we set the flag to one
   options.EVFlag = 1; 
end


%------------------------------------------------------------------------------
function [s,msg,msgobj] = psd_parse(varargin)
%PSD_PARSE Parse trailing inputs for PSD functions.
%  [S,MSG] = PSD_PARSE(varargin) parses the input argument list and
%  returns a structure, S, with fields corresponding to each of the
%  possible PSD options.  If an option is omitted, a default value
%  is returned such that the structure will always return a value
%  for every possible option field.
%
%  MSG is a error message returned from the parse operation.
%  It may be empty, in which case no parsing errors have occurred.
%
%  The structure fields are as follows:
%   S.NFFT    - FFT length (scalar).  If it appears in the argument
%               list, it must be the 1st scalar argument, but not
%               necessarily the 1st argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%
%   S.Fs      - Sample rate in Hz (scalar).  If it appears in the argument
%               list, it must be the 2nd scalar argument, but not
%               necessarily the 2nd argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%
%   S.ConfInt - Confidence interval (scalar).  If it appears in the argument
%   S.nw        list, it must be the 3rd scalar argument, but not
%   (synonyms)  necessarily the 3rd argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%               NOTE: Only pmusic expects S.NW at present, and S.ConfInt and
%               S.NW are simply synonyms of each other (both are set to exactly
%               the same values).
%
%   S.ConfLevel Confidence level property value (scalar) between 0 and 1.
%               If omitted, 'omitted' is returned.  This property does not
%               interfere with S.ConfInt which is used only with PMTM.
%
%   S.DCFlag  - CenterDC.  Returns true iff 'centered' is specified in the argument
%               list.  Recognized argument strings: 'centered'
%
%   S.noverlap  Sample overlap (scalar).  If it appears in the argument list,
%               it must be the 4th scalar argument, but not necessarily the
%               4th argument in the list.  A string option could occur before,
%               instead of, or after this option.  If omitted, 'omitted' is
%               returned.  NOTE: Only pmusic expects S.NOverlap at present.
%
%   S.Range   - X-axis scaling method (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings:
%                       Input             Field Contents
%                       'half'            'one-sided'
%                       'whole'           'two-sided'
%                       'one-sided'       'one-sided'
%                       'two-sided'       'two-sided'
%               NOTE: The hyphens may be omitted or replaced with a single space.
%
%   S.MTMethod -MTM method name (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'adapt', 'eigen', 'unity'.
%
%   S.CorrFlag -Method for pmusic (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'corr'
%
%  s.EVFlag   - Eigenvector method for pmusic (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'EV'


% Initialize return arguments
msg        = '';
msgobj     = [];
omit       = 'omitted';
s.NFFT     = omit;
s.Fs       = omit;
s.ConfInt  = omit;
s.Range    = omit;
s.MTMethod = omit;
s.nw       = omit;
s.noverlap  = omit;
s.CorrFlag = omit;
s.EVFlag   = omit;
s.DCFlag   = omit;
s.ConfLevel = omit;

% look for conflevel pv pair and set flag if user-specified
matchIdx = strcmpi('ConfidenceLevel',varargin);
numMatches = sum(matchIdx(:));
if numMatches>1
   msgobj = message('signal:psdoptions:MultipleValues');
   msg = getString(msgobj);
   return
elseif numMatches == 1
   pIdx = find(matchIdx);
   if pIdx == numel(varargin)
      msgobj = message('signal:psdoptions:MissingConfLevelValue');
      msg = getString(msgobj);
      return
   end
   % obtain the property value.
   confLevel = varargin{pIdx+1};
   if isscalar(confLevel) && isreal(confLevel) && 0<confLevel && confLevel<1
      s.ConfLevel = confLevel;
   else
      msgobj = message('signal:psdoptions:InvalidConfLevelValue');
      msg = getString(msgobj);
      return
   end
   % remove the pv-pair from the argument list
   varargin([pIdx pIdx+1]) = [];    
end

% look for flags

% Define all possible string options
% Lower case, no punctuation:
strOpts = {'half','onesided','whole','twosided', ...
           'adapt','unity','eigen', ...
           'corr', 'ev', ...
           'centered'};

% Check for mutually exclusive options
exclusiveOpts = strOpts([1 3; 2 4; 5 6; 5 7; 6 7; 1 10; 2 10]);
for i=1:size(exclusiveOpts,1)
   if any(strcmpi(exclusiveOpts(i,1),varargin)) && ...
         any(strcmpi(exclusiveOpts(i,2),varargin))
     msgobj = message('signal:psdoptions:ConflictingOptions', ...
                      exclusiveOpts{i,1}, exclusiveOpts{i,2});
     msg = getString(msgobj);
     return
   end
end

% No options passed - return all defaults:
if numel(varargin)<1, return; end

numReals = 0;  % Initialize count of numeric options

for argNum=1:numel(varargin),
   arg = varargin{argNum};
   
   isStr     = ischar(arg);
   isRealVector = isnumeric(arg) & ~issparse(arg) & isreal(arg) & any(size(arg)<=1);
   isValid   = isStr | isRealVector;
   
   if ~isValid,
      msgobj = message('signal:psdoptions:NeedValidOptionsType');
      msg = getString(msgobj);
      return;
   end
   
   if isStr,
      % String option
      %
      % Convert to lowercase and remove special chars:
      arg = RemoveSpecialChar(arg);
      
      % Match option set:
      i = find(strncmp(arg, strOpts, length(arg)));
      if isempty(i) || length(i)>1,
         msgobj = message('signal:psdoptions:UnknStringOption');
         msg = getString(msgobj);
         return
      end
      
      % Determine option set to which this string applies:
      switch i
      case {1,2,3,4},
         field = 'Range';
      case {5,6,7},
         field = 'MTMethod';
      case 8,
         field = 'CorrFlag';
      case 9,
         field = 'EVFlag';
      case 10,
         field = 'DCFlag';
      otherwise,
         error(message('signal:psdoptions:IdxOutOfBound'));
      end
      
      % Verify that no other options from (exclusive) set have
      % already been parsed:
      existingOpt = s.(field);
      if ~strcmp(existingOpt,'omitted'),
         msgobj = message('signal:psdoptions:MultipleValues');
         msg = getString(msgobj);
         return
      end
      
      % Map half to one-sided (index 1 -> index 2)
      % Map whole to two-sided (index 3 -> index 4)
      if i==1 || i==3,
         i=i+1;
      end
      
      % Set full option string into appropriate field:
      s.(field) = strOpts{i};
      
   else
      % Non-string options
      %
      numReals = numReals + 1;
      switch numReals
      case 1, s.NFFT    = arg;
      case 2,
         if length(arg)<=1,
            s.Fs = arg;
         else
            msgobj = message('signal:psdoptions:FsMustBeScalar');
            msg = getString(msgobj);
            return
         end
      case 3,
            s.ConfInt = arg;
            s.nw = arg;  % synonym
            
            % NOTE: Cannot error out for non-scalars, as NW may be a vector!
            %       We cannot disambiguate ConfInt from NW strictly from context.
            %
            %if length(arg)<=1,
            %   s.ConfInt = arg;
            %else
            %   msg = 'The confidence interval must be a scalar.';
            %   return
            %end
      case 4,
         if length(arg)<=1,
            s.noverlap = arg;
         else
            msgobj = message('signal:psdoptions:OverlapMustBeScalar');
            msg = getString(msgobj);
            return
         end
      otherwise
         msgobj = message('signal:psdoptions:TooManyNumericOptions');
         msg = getString(msgobj);
         return
      end
   end
   
end


%--------------------------------------------------------------------------------
function y=RemoveSpecialChar(x)
% RemoveSpecialChar
%   Remove one space of hyphen from 4th position, but
%   only if first 3 chars are 'one' or 'two'

y = lower(x);

% If string is less than 4 chars, do nothing:
if length(y)<=3, return; end

% If first 3 chars are not 'one' or 'two', do nothing:
if ~strncmp(y,'one',3) && ~strncmp(y,'two',3), return; end

% If 4th char is space or hyphen, remove it
if y(4)==' ' || y(4) == '-', y(4)=''; end

% [EOF] psdoptions.m

function p = nextpow2(n)
%NEXTPOW2 Next higher power of 2.
%   NEXTPOW2(N) returns the first P such that 2.^P >= abs(N).  It is
%   often useful for finding the nearest power of two sequence
%   length for FFT operations.
%
%   Class support for input N:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also LOG2, POW2.

%   Copyright 1984-2012 The MathWorks, Inc. 

if ~isinteger(n)
  [f,p] = log2(abs(n));

  % Check for exact powers of 2.
  k = (f == 0.5);
  p(k) = p(k)-1;

  % Check for infinities and NaNs
  k = ~isfinite(f);
  p(k) = f(k);

else % integer case
  p = zeros(size(n),class(n));
  nabs = abs(n);
  x = bitshift(nabs,-1);
  while any(x(:))
    p = p + sign(x);
    x = bitshift(x,-1);
  end
  % Adjust for all non powers of 2
  p = p + max(0,sign(nabs - bitshift(ones(class(n)),p)));
end

function [x,M,isreal_x,y,Ly,win,winName,winParam,noverlap,k,L,options] = ...
    welchparse(x,esttype,varargin)
%WELCHPARSE   Parser for the PWELCH & SPECTROGRAM functions
%
% Outputs:
% X        - First input signal (used when esttype = MS & PSD)
% M        - An integer containing the length of the data to be segmented
% isreal_x - Boolean for input complexity
% Y        - Second input signal (used when esttype = CPSD, TFE, MSCOHERE)
% Ly       - Length of second input signal (used when esttype = CPSD, TFE,
%          MSCOHERE)
% WIN - A scalar or vector containing the length of the window or the
%       window respectively (Note that the length of the window determines
%       the length of the segments)
% WINNAME  - String with the window name.
% WINPARAM - Window parameter.
% NOVERLAP - An integer containing the number of samples to overlap (may
%          be empty)
% K        - Number of segments
% OPTIONS  - A structure with the following fields:
%   OPTIONS.nfft  - number of freq. points at which the psd is estimated
%   OPTIONS.Fs    - sampling freq. if any
%   OPTIONS.range - 'onesided' or 'twosided' psd

%   Copyright 1988-2014 The MathWorks, Inc.

% Parse input arguments.
[x,M,isreal_x,y,Ly,win,winName,winParam,noverlap,opts] = ...
    parse_inputs(x,esttype,varargin{:});

% Obtain the necessary information to segment x and y.
[L,noverlap,win] = segment_info(M,win,noverlap);

% Parse optional args nfft, fs, and spectrumType.
options = welch_options(isreal_x,L,opts{:});

% Compute the number of segments
k = (M-noverlap)./(L-noverlap);

% Uncomment the following line to produce a warning each time the data
% segmentation does not produce an integer number of segments.
%if fix(k) ~= k),
%   warning('signal:welchparse:MustBeInteger','The number of segments is not an integer, truncating data.');
%end

k = fix(k);

%-----------------------------------------------------------------------------------------------
function [x,Lx,isreal_x,y,Ly,win,winName,winParam,noverlap,opts] = ...
    parse_inputs(x,esttype,varargin)
% Parse the inputs to the welch function.

% Assign defaults in case of early return.
y        = [];
Ly       = 0;
is2sig   = false;
win      = [];
winName  = 'User Defined';
winParam = '';
noverlap = [];
opts     = {};

% Determine if one or two signal vectors was specified.
if iscell(x)
    if numel(x) > 1, % Cell array.
        y = x{2};
        is2sig = true;
    end
    x = x{1};
else
    if ~any(strcmpi(esttype,{'psd','power','ms'}))
        error(message('signal:welchparse:NeedCellArray'));
    end
end

if isvector(x)
  x = x(:);
end
Lx = size(x,1);

isreal_x = isreal(x);

% Parse 2nd input signal vector.
if is2sig
    if isvector(y)
      y = y(:);
    end
    isreal_x = isreal(y) && isreal_x;
    Ly = size(y,1);
    if size(x,1) ~= size(y,1)
        error(message('signal:welchparse:MismatchedLength'));
    end
    if size(x,2)~=1 && size(y,2) ~=1 && size(x,2) ~= size(y,2)
        error(message('signal:welchparse:MismatchedNumberOfChannels'));
    end
end

% Parse window and overlap, and cache remaining inputs.
lenargin = length(varargin);
if lenargin >= 1
    win = varargin{1};
    if lenargin >= 2
        noverlap = varargin{2};
        
        % Cache optional args nfft, fs, and spectrumType.
        if lenargin >= 3,  opts = varargin(3:end); end
    end
end

if isempty(win) || isscalar(win)
    winName = 'hamming';
    winParam = 'symmetric';
end

%-----------------------------------------------------------------------------------------------
function [L,noverlap,win] = segment_info(M,win,noverlap)
%SEGMENT_INFO   Determine the information necessary to segment the input data.
%
%   Inputs:
%      M        - An integer containing the length of the data to be segmented
%      WIN      - A scalar or vector containing the length of the window or the window respectively
%                 (Note that the length of the window determines the length of the segments)
%      NOVERLAP - An integer containing the number of samples to overlap (may be empty)
%
%   Outputs:
%      L        - An integer containing the length of the segments
%      NOVERLAP - An integer containing the number of samples to overlap
%      WIN      - A vector containing the window to be applied to each section
%
%
%   The key to this function is the following equation:
%
%      K = (M-NOVERLAP)/(L-NOVERLAP)
%
%   where
%
%      K        - Number of segments
%      M        - Length of the input data X
%      NOVERLAP - Desired overlap
%      L        - Length of the segments
%
%   The segmentation of X is based on the fact that we always know M and two of the set
%   {K,NOVERLAP,L}, hence determining the unknown quantity is trivial from the above
%   formula.

% Initialize outputs
L = [];

% Check that noverlap is a scalar
if any(size(noverlap) > 1)
    error(message('signal:welchparse:invalidNoverlap'));
end

if isempty(win)
    % Use the closest to 8 sections, determine their length
    if isempty(noverlap)
        % Use 50% overlap
        L = fix(M./4.5);
        noverlap = fix(0.5.*L);
    else
        L = fix((M+7.*noverlap)./8);
    end
    % Use a default window
    win = hamming(L);

else
    % Determine the window and its length (equal to the length of the segments)
    if ~any(size(win) <= 1) || ischar(win)
        error(message('signal:welchparse:MustBeScalarOrVector', 'WINDOW'));
    elseif length(win) > 1
        % WIN is a vector
        L = length(win);
    elseif length(win) == 1
        L = win;
        win = hamming(win);
    end
    if isempty(noverlap)
        % Use 50% overlap
        noverlap = fix(0.5.*L);
    end
end

% Do some argument validation
if L > M
    error(message('signal:welchparse:invalidSegmentLength'));
end

if noverlap >= L
    error(message('signal:welchparse:NoverlapTooBig'));
end

%------------------------------------------------------------------------------
function options = welch_options(isreal_x,N,varargin)
%WELCH_OPTIONS   Parse the optional inputs to the PWELCH function.
%   WELCH_OPTIONS returns a structure, OPTIONS, with following fields:
%
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided' psd

% Generate defaults
options.nfft = max(256,2^nextpow2(N));
options.Fs = []; % Work in rad/sample

% Determine if frequency vector specified
freqVecSpec = false;
if (~isempty(varargin) && length(varargin{1}) > 1)
    freqVecSpec = true;
end

if isreal_x && ~freqVecSpec
    options.range = 'onesided';
else
    options.range = 'twosided';
end

if any(strcmp(varargin, 'whole'))
    warning(message('signal:welchparse:InvalidRange', '''whole''', '''twosided'''));
elseif any(strcmp(varargin, 'half'))
    warning(message('signal:welchparse:InvalidRange', '''half''', '''onesided'''));
end

[options,msg,msgobj] = psdoptions(isreal_x,options,varargin{:});
if ~isempty(msg), error(msgobj); end;



% [EOF]

function varargout = welch(x,esttype,varargin)
%WELCH Welch spectral estimation method.
%   [Pxx,F] = WELCH(X,WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,ESTTYPE)
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'psd')
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'ms')
%   [Pxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'cpsd')
%   [Txy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'tfe')
%   [Cxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'mscohere')
%   [Pxx,F,Pxxc] = WELCH(...)
%   [Pxx,F,Pxxc] = WELCH(...,'ConfidenceLevel',P)
%
%   Inputs:
%      see "help pwelch" for complete description of all input arguments.
%      ESTTYPE - is a string specifying the type of estimate to return, the
%                choices are: psd, cpsd, tfe, and mscohere.
%
%   Outputs:
%      Depends on the input string ESTTYPE:
%      Pxx - Power Spectral Density (PSD) estimate, or
%      MS  - Mean-square spectrum, or
%      Pxy - Cross Power Spectral Density (CPSD) estimate, or
%      Txy - Transfer Function Estimate (TFE), or
%      Cxy - Magnitude Squared Coherence.
%      F   - frequency vector, in Hz if Fs is specified, otherwise it has
%            units of rad/sample

%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and
%         Modeling, John Wiley & Sons, 1996.

narginchk(2,10);
nargoutchk(0,3);

% Parse input arguments.
[x,~,~,y,~,win,winName,winParam,noverlap,k,L,options] = ...
    welchparse(x,esttype,varargin{:});
% Cast to enforce precision rules
options.nfft = sigcasttofloat(options.nfft,'double',...
  'WELCH','NFFT','allownumeric');
noverlap = sigcasttofloat(noverlap,'double','WELCH',...
  'NOVERLAP','allownumeric');
options.Fs = sigcasttofloat(options.Fs,'double','WELCH',...
  'Fs','allownumeric');
k = double(k);

if any([sigcheckfloattype(x,'single')...
    sigcheckfloattype(y,'single'),...
    isa(win,'single')])
  x = single(x);
  y = single(y);
  win = single(win);
end

% Frequency vector was specified, return and plot two-sided PSD
freqVectorSpecified = false; nrow = 1;
if length(options.nfft) > 1,
    freqVectorSpecified = true;
    [~,nrow] = size(options.nfft);
end

% Compute the periodogram power spectrum of each segment and average always
% compute the whole power spectrum, we force Fs = 1 to get a PS not a PSD.

% Initialize
if freqVectorSpecified,
    nFreqs = length(options.nfft);
else
    nFreqs = options.nfft;
end
% Cast to enforce precision rules
Sxx = zeros(nFreqs,size(x,2),class(x)); %#ok<*ZEROLIKE>

LminusOverlap = L-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+L-1;
switch esttype,
    case {'ms','power','psd'}
        for i = 1:k
            [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i),:),win,...
                options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
        end
        
    case 'cpsd'
        numChan = max(size(x,2),size(y,2));
        Sxx = zeros(nFreqs,numChan,class(x)); % Initialize
        for i = 1:k
            [Sxxk,w] =  computeperiodogram({x(xStart(i):xEnd(i),:),...
                y(xStart(i):xEnd(i),:)},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
        end
        
    case 'tfe'
        numChan = max(size(x,2),size(y,2));
        Sxy = zeros(nFreqs,numChan,class(x));
        for i = 1:k
            Sxxk = computeperiodogram(x(xStart(i):xEnd(i),:),...
                win,options.nfft,esttype,options.Fs);
            [Syxk,w] = computeperiodogram({y(xStart(i):xEnd(i),:),...
                x(xStart(i):xEnd(i),:)},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
            Sxy  = Sxy + Syxk;
        end
        
    case 'mscohere'
        % Note: (Sxy1+Sxy2)/(Sxx1+Sxx2) != (Sxy1/Sxy2) + (Sxx1/Sxx2)
        % ie, we can't push the computation of Cxy into computeperiodogram.
        numChan = max(size(x,2),size(y,2));
        Sxy = zeros(nFreqs,numChan,class(x));
        Syy = zeros(nFreqs,size(y,2),class(x));
        
        for i = 1:k
            Sxxk = computeperiodogram(x(xStart(i):xEnd(i),:),...
                win,options.nfft,esttype,options.Fs);
            Syyk = computeperiodogram(y(xStart(i):xEnd(i),:),...
                win,options.nfft,esttype,options.Fs);
            [Sxyk,w] = computeperiodogram({x(xStart(i):xEnd(i),:),...
                y(xStart(i):xEnd(i),:)},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
            Syy  = Syy + Syyk;
            Sxy  = Sxy + Sxyk;
        end
end
Sxx = Sxx./k; % Average the sum of the periodograms

if any(strcmpi(esttype,{'tfe','mscohere'})),
    Sxy = Sxy./k; % Average the sum of the periodograms
    
    if strcmpi(esttype,'mscohere'),
        Syy = Syy./k; % Average the sum of the periodograms
    end
end

% Generate the freq vector directly in Hz to avoid roundoff errors due to
% conversions later.
if ~freqVectorSpecified
    w = psdfreqvec('npts',options.nfft, 'Fs',options.Fs);
else
    if strcmpi(options.range,'onesided')
        warning(message('signal:welch:InconsistentRangeOption'));
    end
    options.range = 'twosided';
end


% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, corresponding freq vector and freq units.
[Pxx,w,units] = computepsd(Sxx,w,options.range,options.nfft,options.Fs,esttype);

if any(strcmpi(esttype,{'tfe','mscohere'}))
    % Cross PSD.  The frequency vector and xunits are not used.
    Pxy = computepsd(Sxy,w,options.range,options.nfft,options.Fs,esttype);
    
    % Transfer function estimate.
    if strcmpi(esttype,'tfe')
        Pxx = bsxfun(@rdivide, Pxy, Pxx);
    end
    
    % Magnitude Square Coherence estimate.
    if strcmpi(esttype,'mscohere')
        % Auto PSD for 2nd input vector. The freq vector & xunits are not
        % used.
        Pyy = computepsd(Syy,w,options.range,options.nfft,options.Fs,esttype);
        Pxx = (abs(Pxy).^2)./bsxfun(@times,Pxx,Pyy); % Cxy
    end
end

% compute confidence intervals if needed.
if ~strcmp(options.conflevel,'omitted')
    Pxxc = confInterval(options.conflevel, Pxx, x, w, options.Fs, k);
elseif nargout>2
    Pxxc = confInterval(0.95, Pxx, x, w, options.Fs, k);
else
    Pxxc = [];
end

if nargout==0
    w = {w};
    if strcmpi(units,'Hz'), w = [w,{'Fs',options.Fs}];  end
    % Create a spectrum object to store in the Data object's metadata.
    percOverlap = (noverlap/L)*100;
    hspec = spectrum.welch({winName,winParam},L,percOverlap);
    
    switch lower(esttype)
        case 'tfe'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.freqz(Pxx,w{:},'SpectrumRange',range);
        case 'mscohere'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.magnitude(Pxx,w{:},'SpectrumRange',range);
        case 'cpsd'
            h = dspdata.cpsd(Pxx,w{:},'SpectrumType',options.range);
        case {'ms','power'}
            h = dspdata.msspectrum(Pxx,w{:},'SpectrumType',options.range);
        otherwise
            h = dspdata.psd(Pxx,w{:},'SpectrumType',options.range);
    end
    h.Metadata.setsourcespectrum(hspec);
    
    % plot the confidence levels if conflevel is specified.
    if ~isempty(Pxxc)
        h.ConfLevel = options.conflevel;
        h.ConfInterval = Pxxc;
    end
    % center dc component if specified
    if options.centerdc
        centerdc(h);
    end
    plot(h);
    if strcmp(esttype,'power')
        title(getString(message('signal:welch:WelchPowerSpectrumEstimate')));
    end
else
    if options.centerdc
        [Pxx, w, Pxxc] = psdcenterdc(Pxx, w, Pxxc, options);
    end
    
    % If the input is a vector and a row frequency vector was specified,
    % return output as a row vector for backwards compatibility.
    if nrow > 1 && isvector(x)
        Pxx = Pxx.'; w = w.';
    end
    
    % Cast to enforce precision rules   
    % Only cast if output is requested, otherwise, plot using double
    % precision frequency vector.
    if isa(Pxx,'single')
      w = single(w);
    end
    
    if isempty(Pxxc)
        varargout = {Pxx,w}; % Pxx=PSD, MEANSQUARE, CPSD, or TFE
    else
        varargout = {Pxx,w,Pxxc};
    end       
end

function Pxxc = confInterval(CL, Pxx, x, w, fs, k)
%   Reference: D.G. Manolakis, V.K. Ingle and S.M. Kagon,
%   Statistical and Adaptive Signal Processing,
%   McGraw-Hill, 2000, Chapter 5
k = fix(k);
c = privatechi2conf(CL,k);

% Cast to enforce precision rules
Pxx = double(Pxx);

PxxcLower = Pxx*c(1);
PxxcUpper = Pxx*c(2);
Pxxc = reshape([PxxcLower; PxxcUpper],size(Pxx,1),2*size(Pxx,2));

% DC and Nyquist bins have only one degree of freedom for real signals
if isreal(x)
    realConf = chi2conf(CL,k/2);
    Pxxc(w == 0,1:2:end) = Pxx(w == 0,:) * realConf(1);
    Pxxc(w == 0,2:2:end) = Pxx(w == 0,:) * realConf(2);
    if isempty(fs)
        Pxxc(w == pi,1:2:end) = Pxx(w == pi,:) * realConf(1);
        Pxxc(w == pi,2:2:end) = Pxx(w == pi,:) * realConf(2);
    else
        Pxxc(w == fs/2,1:2:end) = Pxx(w == fs/2,:) * realConf(1);
        Pxxc(w == fs/2,2:2:end) = Pxx(w == fs/2,:) * realConf(2);
    end
end

% [EOF]


