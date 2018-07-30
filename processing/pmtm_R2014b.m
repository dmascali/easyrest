function varargout = pmtm_R2014b(x,varargin)
%PMTM   Power Spectral Density (PSD) estimate via the Thomson multitaper 
%   method (MTM).
%   Pxx = PMTM(X) returns the Power Spectral Density (PSD) estimate, Pxx,
%   of a discrete-time signal X. When X is a vector, it is converted to a
%   column vector and treated as a single channel.  When X is a matrix, the
%   PSD is computed independently for each column and stored in the
%   corresponding column of Pxx. Pxx is the distribution of power per unit
%   frequency. The frequency is expressed in units of radians/sample.
%
%   For real signals, PMTM returns the one-sided PSD by default; for 
%   complex signals, it returns the two-sided PSD.  Note that a one-sided 
%   PSD contains the total power of the input signal.
%
%   Pxx = PMTM(X,NW) specifies NW as the "time-bandwidth product" for the
%   discrete prolate spheroidal sequences (or Slepian sequences) used as 
%   data windows.  Typical choices for NW are 2, 5/2, 3, 7/2, or 4. If
%   empty or omitted, NW defaults to 4. By default, PMTM drops the last
%   taper because its corresponding eigenvalue is significantly smaller
%   than 1. Therefore, The number of tapers used to form Pxx is 2*NW-1.
%   
%   Pxx = PMTM(X,NW,NFFT) specifies the FFT length used to calculate the
%   PSD estimates.  For real X, Pxx has (NFFT/2+1) rows if NFFT is even,
%   and (NFFT+1)/2 rows if NFFT is odd.  For complex X, Pxx always has
%   length NFFT.  If NFFT is specified as empty, NFFT is set to either
%   256 or the next power of 2 greater than the length of X, whichever is
%   larger.
%
%   [Pxx,W] = PMTM(X,NW,NFFT) returns the vector of normalized angular 
%   frequencies, W, at which the PSD is estimated.  W has units of 
%   radians/sample.  For real signals, W spans the interval [0,Pi] when
%   NFFT is even and [0,Pi) when NFFT is odd.  For complex signals, W 
%   always spans the interval [0,2*Pi).
%
%   [Pxx,W] = PMTM(X,NW,W) uses the Goertzel algorithm to compute the
%   two-sided PSD at the normalized angular frequencies contained in vector
%   W. W must have at least two elements. The specified frequencies in W
%   are rounded to the nearest DFT bin commensurate with the signal's
%   resolution.
%
%   [Pxx,F] = PMTM(X,NW,NFFT,Fs) returns a PSD computed as a function of
%   physical frequency.  Fs is the sampling frequency specified in hertz.
%   If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which the PSD is
%   estimated.  For real signals, F spans the interval [0,Fs/2] when NFFT
%   is even and [0,Fs/2) when NFFT is odd.  For complex signals, F always
%   spans the interval [0,Fs).
%
%   [Pxx,F] = PMTM(X,NW,F,Fs) where F is a vector of 
%   frequencies in Hz (with 2 or more elements) computes the PSD at 
%   those frequencies using the Goertzel algorithm. In this case a two
%   sided PSD is returned. The specified frequencies in F are rounded to 
%   the nearest DFT bin commensurate with the signal's resolution. 
%
%   [Pxx,F] = PMTM(...,Fs,method) uses the algorithm specified in method 
%   for combining the individual spectral estimates:
%      'adapt'  - Thomson's adaptive non-linear combination (default).
%      'unity'  - linear combination with unity weights.
%      'eigen'  - linear combination with eigenvalue weights.
%
%   [Pxx,F] = PMTM(X,E,V,NFFT,Fs,method) is the PSD estimate, and frequency
%   vector from the data tapers in E and their concentrations V.  Type HELP
%   DPSS for a description of the matrix E and the vector V. By default,
%   PMTM drops the last eigenvector because its corresponding eigenvalue is
%   significantly smaller than 1.
%
%   [Pxx,F] = PMTM(X,DPSS_PARAMS,NFFT,Fs,method) uses the cell 
%   array DPSS_PARAMS containing the input arguments to DPSS (listed in
%   order, but excluding the first argument) to compute the data tapers. 
%   For example, PMTM(x,{3.5,'trace'},512,1000) calculates the prolate 
%   spheroidal sequences for NW=3.5, NFFT=512, and Fs=1000, and displays
%   the method that DPSS uses for this calculation. Type HELP DPSS for 
%   other options.
%
%   [Pxx,F] = PMTM(...,'DropLastTaper',DROPFLAG) specifies whether PMTM
%   should drop the last taper/eigenvector during the calculation. DROPFLAG
%   can be one of the following values: [ {true} | false ].
%       true  - the last taper/eigenvector is dropped 
%       false - the last taper/eigenvector is preserved
%
%   [Pxx,F,Pxxc] = PMTM(...,'ConfidenceLevel',P) returns the P*100%
%   confidence interval for Pxx, where P is a scalar between 0 and 1. The
%   default value for P is .95.  Confidence intervals are computed using a
%   chi-squared approach. Pxxc has twice as many columns as Pxx.
%   Odd-numbered columns contain the lower bounds of the confidence
%   intervals; even-numbered columns contain the upper bounds.  Thus,
%   Pxxc(M,2*N-1) is the lower bound and Pxxc(M,2*N) is the upper bound
%   corresponding to the estimate Pxx(M,N).
%
%   [...] = PMTM(X,...,FREQRANGE) returns the PSD over the specified range
%   of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided PSD of a real input signal X.
%         If NFFT is even, Pxx has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, Pxx has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided PSD for either real or complex
%         input X.  Pxx has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided PSD for either real or
%         complex X.  Pxx has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd length NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after the second input argument, unless E and V are specified, in
%      which case FREQRANGE may be placed in any position after the third
%      input argument.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   PMTM(...) with no output arguments plots the PSD estimate (in decibels
%   per unit frequency) in the current figure window.
%
%   EXAMPLE:
%      Fs = 1000;   t = 0:1/Fs:.3;  
%      x = cos(2*pi*t*200)+randn(size(t)); % A cosine of 200Hz plus noise
%      pmtm(x,3.5,[],Fs);                  % Uses the default NFFT.
%
%   See also DPSS, PWELCH, PERIODOGRAM, PMUSIC, PBURG, PYULEAR, PCOV,
%   PMCOV, PEIG.

%   References: 
%     [1] Thomson, D.J."Spectrum estimation and harmonic analysis."
%         In Proceedings of the IEEE. Vol. 10 (1982). pp. 1055-1096.
%     [2] Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
%         Applications", Cambridge University Press, 1993, pp. 368-370. 

%   Copyright 1988-2014 The MathWorks, Inc.
%      

narginchk(1,13);
len = length(varargin);

% convert vectors to column vectors for backwards compatibility
if isvector(x)
  x = x(:);
end

if numel(size(x))>2
  error(message('signal:pmtm:MustBeVectorOr2DMatrix'));
end

tf = strcmpi('droplasttaper',varargin);
indx = find(tf==1);
if (~isempty(indx) && ~islogical(varargin{indx+1}))
    error(message('signal:pmtm:MustBeLogical', 'Droplasttaper'));
end

% If the 'droplasttaper' pv-pair is used, move it to the end of varargin
if (~isempty(indx) && (indx+1 ~= len))
    dummy = varargin(1:indx-1);
    dummy(indx:len-2) = varargin(indx+2:len);
    dummy(len-1:len) = varargin(indx:indx+1);
    varargin = dummy;
end

% Parse the inputs, set up default values, and return any error messages.
params = parseinputs(x,varargin{:});

% Cast to enforce Precision Rules
% Already checked for invalid character inputs (NW, NFFT,Fs) in 'parseinputs->psdoptions'
params.nfft = double(params.nfft);
params.Fs = double(params.Fs);

% Compute the two-sided power spectrum via MTM.
[S,k,w] = mtm_spectrum(x,params);

% Generate the freq vector in correct units to avoid roundoff errors due to
% converting frequency units later.
nfft = params.nfft;
[~,ncol] = size(nfft);

% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, compute the corresponding freq vector & freq units.
[Pxx,f,units] = computepsd(S,w,params.range,params.nfft,params.Fs,'psd');  

% Calculate confidence limits ONLY when needed, since it can take a while.
% Enforce double precision arithmetic on the calculations. chi2conf already
% enforces double precision arithmetic.
if ~strcmp(params.conflevel,'omitted') && (nargout==0 || nargout>2)
  % 'ConfidenceLevel' pv-pair specified
  Nchan = size(x,2);
  c = chi2conf(params.conflevel,k);
  Pxxc(:,2:2:2*Nchan) = double(Pxx)*c(2);
  Pxxc(:,1:2:2*Nchan) = double(Pxx)*c(1);
elseif (~strcmp(params.ConfInt,'omitted') && nargout==0) || nargout>2
  % using legacy syntax
  if ~strcmp(params.ConfInt,'omitted')
    c = chi2conf(params.ConfInt,k);
  else % use default value of .95
    c = chi2conf(.95,k);
  end      
  Nchan = size(x,2);
  Pxxc(:,2:2:2*Nchan) = double(Pxx)*c(2);
  Pxxc(:,1:2:2*Nchan) = double(Pxx)*c(1);
else
  Pxxc = [];
end

% Output
if nargout==0
   % If no output arguments are specified plot the PSD w/ conf intervals.
   f = {f};
   if strcmpi(units,'Hz'), f = [f {'Fs',params.Fs}]; end
   hpsd = dspdata.psd([Pxx Pxxc],f{:},'SpectrumType',params.range);

   % Create a spectrum object to store in the PSD object's metadata.
   hspec = spectrum.mtm(params.E,params.V,params.MTMethod);
   hpsd.Metadata.setsourcespectrum(hspec);
      
   if params.centerdc
     centerdc(hpsd);
   end

   hp = plot(hpsd);
   if 3*size(x,2)==numel(hp)
     nChan = size(x,2);
     color = get(hp,'Color');
     for i=1:nChan
       set(hp(nChan+2*i-1), 'Color',color{i}, 'LineStyle','-.');
       set(hp(nChan+2*i), 'Color',color{i}, 'LineStyle','-.');
     end
   end     
   return
end

% center dc if needed
if params.centerdc
   [Pxx, f, Pxxc] = psdcenterdc(Pxx, f, Pxxc, params);
end

% If the input is a vector and a row frequency vector was specified,
% return output as a row vector for backwards compatibility.
if ncol > 1 && nargout > 0 && isvector(x)
    Pxx = Pxx.';
    f = f.';
    % preserve (undocumented) behavior with legacy syntax.
    if strcmp(params.conflevel,'omitted') && nargout >= 3
        Pxxc = Pxxc.';
    end
end

if isa(Pxx,'single')
  % Cast to enforce precision rules.
  f = single(f);
end  

if nargout==1
   varargout = {Pxx};
elseif nargout==2
   varargout = {Pxx,f};
elseif nargout==3
   if ~strcmp(params.conflevel,'omitted')
      % use preferred output order
      varargout = {Pxx,f,Pxxc};
   else
      % use legacy output order
      varargout = {Pxx,Pxxc,f};
   end
end

%----------------------------------------------------------------------
function [S,k,w] = mtm_spectrum(x,params)
%MTM_SPECTRUM Compute the power spectrum via MTM.
%
% Inputs:
%   x      - Input data vector.
%   params - Structure containing pmtm's input parameter list, except for
%            the input data sequence, x; it contains the following fields:
%      nfft     - Number of frequency points to evaluate the PSD at; 
%                 the default is max(256,2^nextpow2(N)).
%      Fs       - The sampling frequency; default is 1.
%      range    - default is 'onesided' or real signals and 'twosided' for 
%               - complex signals.
%      ConfInt  - Confidence interval; default is .95.
%      MTMethod - Algorithm used in MTM; default is 'adapt'.
%      E        - Matrix containing the discrete prolate spheroidal 
%                 sequences (dpss).
%      V        - Vector containing the concentration of the dpss.
%      NW       - Time-bandwidth product; default is 4.
%
% Outputs:
%   S      - Power spectrum computed via MTM.
%   k      - Number of sequences used to form Pxx.
%   w      - Frequency vector for which DFT is calculated

% Extract some parameters from the input structure for convenience.
nfft = params.nfft;
E  = params.E;
V  = params.V;
Fs = params.Fs;

if isempty(Fs)
  Fs = 2*pi;
end

N = size(x,1);
Nchan = size(x,2);
k = length(V);

if length(nfft) > 1, 
    isfreqVector = true;     
    nfft_mod = length(nfft);
else 
    isfreqVector = false;
    nfft_mod = nfft;
end

if isa(x,'single') || isa(E,'single')
  S = zeros(nfft_mod, Nchan, 'single');
else
  S = zeros(nfft_mod, Nchan);
end

for chan=1:Nchan
    % Compute the windowed DFTs.
    if (~isfreqVector && N<=nfft) || isfreqVector 

        % Compute DFT using FFT or Goertzel
        [Xx,w] = computeDFT(bsxfun(@times,E(:,1:k),x(:,chan)),nfft(:),Fs);    
        Sk = abs(Xx).^2;

    else % Wrap the data modulo nfft if N > nfft. Note we cannot use datawrap 
        % and FFT because datawrap does not support matrices
        % use CZT to compute DFT on nfft evenly spaced samples around the
        % unit circle:
        Sk = abs(czt(bsxfun(@times,E(:,1:k),x(:,chan)),nfft(:))).^2;
        w = psdfreqvec('npts',nfft,'Fs',Fs);
    end

    % Compute the MTM spectral estimates, compute the whole spectrum 0:nfft.
    switch params.MTMethod,

    case 'adapt'
       % Set up the iteration to determine the adaptive weights: 

       sig2=x(:,chan)'*x(:,chan)/N;              % Power
       Schan=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate   
       S1=zeros(nfft_mod,1);  

       % The algorithm converges so fast that results are
       % usually 'indistinguishable' after about three iterations.

       % This version uses the equations from [2] (P&W pp 368-370).

       % Set tolerance for acceptance of spectral estimate:
       tol=.0005*sig2/nfft_mod;
       a=bsxfun(@times,sig2,(1-V));

       % Do the iteration:
       while sum(abs(Schan-S1)/nfft_mod)>tol
          % calculate weights
          b=(Schan*ones(1,k))./(Schan*V'+ones(nfft_mod,1)*a'); 
          % calculate new spectral estimate
          wk=(b.^2).*(ones(nfft_mod,1)*V');
          S1=sum(wk'.*Sk')./ sum(wk,2)';
          S1=S1';
          Stemp=S1; S1=Schan; Schan=Stemp;  % swap S and S1
       end
    case {'unity','eigen'}
       % Compute the averaged estimate: simple arithmetic averaging is used. 
       % The Sk can also be weighted by the eigenvalues, as in Park et al. 
       % Eqn. 9.; note that the eqn. apparently has a typo; as the weights
       % should be V and not 1/V.
       if strcmp(params.MTMethod,'eigen')
          wt = V(:);    % Park estimate
       else
          wt = ones(k,1);
       end
       Schan = Sk*wt/k;
    end
    S(:,chan) = Schan;
end

%----------------------------------------------------------------------
function params = parseinputs(x,varargin)
%PARSEINPUTS Parse the inputs passed to pmtm.m and return a structure
%            containing all the parameters passed to PMTM set to either
%            default values or user defined values.
%
% Inputs:
%   x        - Input data vector.
%   varargin - Input parameter list passed to pmtm, except for x.
%
% Outputs:
%   params   - Structure containing pmtm's input parameter list, except for
%              the input data sequence, x; it contains the following fields:
%      nfft     - Number of frequency points to evaluate the PSD at; 
%                 the default is max(256,2^nextpow2(N)).
%      Fs       - The sampling frequency; default is .
%      range    - default is 'onesided' or real signals and 'twosided' for 
%               - complex signals.
%      conflevel- Confidence level (preferred syntax)
%      ConfInt  - Confidence interval; default is .95. (legacy syntax)
%      MTMethod - Algorithm used in MTM; default is 'adapt'.
%      E        - Matrix containing the discrete prolate spheroidal 
%                 sequences.
%      V        - Vector containing the concentration of the dpss.
%      NW       - Time-bandwidth product; default is 4.
%
%   err_msg  - String containing an error message if an error occurred.

if any(strcmp(varargin, 'whole'))
    warning(message('signal:pmtm:invalidRange', 'whole', 'twosided'));
elseif any(strcmp(varargin, 'half'))
    warning(message('signal:pmtm:invalidRange', 'half', 'onesided'));
end

% Set default parameter values.
N = size(x,1);
params  = [];

% Parse the input arguments up to NFFT (exclusive). 
% If E and V are not specified, calculate them.
[E,V,NW,indx,nfft_temp,varargin] = getEV(N,varargin{:});

% Cast to enforce Precision Rules
if (any([sigcheckfloattype(x,'single','pmtm','X') ...
    sigcheckfloattype(E,'single','pmtm','E') ...
    sigcheckfloattype(V,'single','pmtm','V')]))
%      sigcheckfloattype(E,'single','pmtm','E') ...
%      sigcheckfloattype(V,'single','pmtm','V')]))

  x = single(x);
  E = single(E);
  V = single(V);
end
NW = double(NW);

if isreal(x) && (length(nfft_temp) <= 1), 
   range = 'onesided';
else
   range = 'twosided'; 
end

% NOTE: The psdoptions function REQUIRES a structure with the following 
%       fields.  Any changes to the structure, such as adding/removing 
%       fields, should be done after the call to psdoptions.
params.nfft    = max(256,2^nextpow2(N));
params.Fs      = [];
params.range   = range;
params.centerdc = false;
params.conflevel = 'omitted';

params.ConfInt = 'omitted';
params.MTMethod= 'adapt';

% Call psdoptions to handle the remaining input arg list starting with NFFT.
% Overwrite default options with user specified options (if specified).
if indx <= numel(varargin)
  % Invalid character inputs for NW, NFFT, W, E,V and Fs is checked here
   [params,err_msg,err_msgobj] = psdoptions(isreal(x),params,varargin{indx:end});
   if err_msg, error(err_msgobj), end;     
   
   if ~strcmp(params.conflevel,'omitted') && ~strcmp(params.ConfInt, 'omitted')
       % user specified too many scalar inputs in conjunction with 'ConfidenceLevel'
       error(message('signal:pmtm:TooManyScalarNumericInputs'));
   end
   
   if length(params.nfft) > 1,
       if strcmpi(params.range,'onesided')
           warning(message('signal:pmtm:InconsistentRangeOption'));
       end
       params.range = 'twosided';
   end
end

% Add remaining fields to the return structure.
params.E  = E;
params.V  = V;
params.NW = NW;

%----------------------------------------------------------------------
function [E,V,NW,indx,nfft_temp,varargin] = getEV(N,varargin)
% GETEV  Parse the input arguments up to, but not including, Nfft and 
%        calculate E and V if not specified.
%
% Inputs:
%   N        - Length of the input data sequence, x.
%   varargin - Input parameter list passed to pmtm, except for x.
%
% Outputs:
%   E        - Matrix containing the discrete prolate spheroidal 
%              sequences (dpss).
%   V        - Vector containing the concentration of the dpss.
%   NW       - Time-bandwidth product; default is 4.
%   indx     - Index indicating starting location of options in pmtm's 
%              input argument list.
%   nfft_temp - NFFT or Frequency vector specified. Empty if not specified 

% Define defaults & initialize output variables (in case of early return).
NW      = 4;
indx    = 1;  % Index where the options begin in the input arg list
nfft_temp = [];

tf = strcmpi('droplasttaper',varargin);
loc = find(tf==1);
if ~isempty(loc)
    dlt = varargin{loc+1};     % droplasttaper
    varargin = varargin(1:loc-1);
else
    dlt = true;              % default value
end

% The 2nd input arg to pmtm can be a
%    1. (X,NW,...)            scalar
%    2. (X,E,V,...)           matrix E, hence, 3rd input must be a vector (V) 
%    3. (X,{dpss_params},...) cell containing the input argument list to dpss 
%    4. a string specifying a psdoption
if ~isempty(varargin) && ~ischar(varargin{1})
   if ~isempty(varargin{1})
       NW = varargin{1};
   end
   indx = 2;
   if iscell(NW),           % NW is a cell array => dpss_params
      if (NW{1}<1.25 && dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '1.25', 'Droplasttaper', 'true'));
      end 
      if (NW{1}<0.75 && ~dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '0.75', 'Droplasttaper', 'false'));          
      end       
      [E,V] = dpss(N,NW{:}); 
      numvec = length(V);
      if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
      else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
      end
      NW = NW{1};
      if nargin > 2, nfft_temp = findNFFT(varargin{2:end}); end
   elseif length(NW)>1,     % NW is the matrix E (==>V must be specified)
      E = NW;
      if length(varargin)<2,
         error(message('signal:pmtm:MustProvideVWithE', 'V', 'E'));
      else
         V = varargin{2};
         if nargin > 3, nfft_temp = findNFFT(varargin{3:end}); end
      end
      numvec = length(V);
      if size(E,2)~=numvec
         error(message('signal:pmtm:MismatchEV', 'E', 'V'));
      end     
      NW = size(E,2)/2;  
      indx = 3; % Update index of beginning of options in the input arg list      
      if dlt
          if (numvec < 3)
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
          else
              E = E(:,1:numvec-1);
              V = V(1:numvec-1);
          end
      else
          if(numvec < 2)
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
          end
      end
   else                      % NW is a scalar
       if (NW<1.25 && dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '1.25', 'Droplasttaper', 'true'));
       end
       if (NW<0.75 && ~dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '0.75', 'Droplasttaper', 'false'));
       end
       % Get the dpss, one way or another:
       [E,V] = dpss(N,NW);
       numvec = length(V);
       if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
       else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
       end
       if nargin > 2, nfft_temp = findNFFT(varargin{2:end}); end
   end
else
   % Get the dpss, one way or another:
   [E,V] = dpss(N,NW);
   numvec = length(V);
      if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
      else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
      end
   nfft_temp = [];
end


%------------------------------------------------------------------
function nfftTemp = findNFFT(varargin)
% FINDNFFT Finds the specified NFFT or frequency vector from the optional
% arguments passed

nfftTemp = [];
for cnt = 1:length(varargin)
    if isnumeric(varargin{cnt}), 
        nfftTemp = varargin{cnt};
        break;
    end
end

% [EOF] pmtm.m

function [E,V]=dpss(N, NW, varargin)
%DPSS   Discrete prolate spheroidal sequences (Slepian sequences).
%   [E,V] = DPSS(N,NW) are the first 2*NW discrete prolate spheroidal sequences
%   (DPSSs, or Slepian sequences) of length N (in the columns of E) and 
%   their corresponding concentrations (in vector V) in the frequency band 
%   |w|<=(2*pi*W) (where  W = NW/N is the half-bandwidth and w is in 
%   radians/sample).  E(:,1) is the length N signal most concentrated in the 
%   frequency band |w|<=(2*pi*W) radians/sample, E(:,2) is the signal 
%   orthogonal to E(:,1) which is most concentrated in this band, E(:,3) is the
%   signal orthogonal to both E(:,1) and E(:,2) which is most concentrated in 
%   this band, etc.  
%
%   For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 
%   7/2, or 4.
%
%   [E,V] = DPSS(N,NW,K) are the K most band-limited discrete prolate spheroidal
%   sequences.  [E,V] = DPSS(N,NW,[K1 K2]) returns the K1-th through the 
%   K2-th sequences.
%
%   [E,V] = DPSS(N,NW,'spline') uses spline interpolation to compute the DPSSs 
%   from existing DPSSs in the DPSS database with length closest to N.
%   [E,V] = DPSS(N,NW,'spline',Ni) interpolates from existing length Ni DPSSs.
%   DPSS(N,NW,'linear') and DPSS(N,NW,'linear',Ni) use linear interpolation, 
%   which is much faster but less accurate than spline interpolation.  
%   'linear' requires Ni > N. [E,V] = DPSS(N,NW,'calc') uses the direct 
%   algorithm (default).
%
%   Use a trailing 'trace' argument to find out which method DPSS uses, e.g.,
%   DPSS(...,'trace').
%
%   % Example:
%   %   Construct a set of Slepian sequences. 
%
%   seq_length = 512; 
%   time_halfbandwidth = 2.5;
%   num_seq = 2*(2.5)-1;
%   %Obtain DPSSs
%   [dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);
%   % Plot the Slepian sequences
%   plot(dps_seq);
%   title('Slepian Sequences N=512, NW=2.5');
%   axis([0 512 -0.15 0.15]);
%   legend('1st','2nd','3rd','4th');
%
%   See also PMTM, DPSSLOAD, DPSSDIR, DPSSSAVE, DPSSCLEAR.

%   References: 
%     [1] Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
%         Applications", Cambridge University Press, 1993. 

%   Author: Eric Breitenberger, 10/3/95
%   Copyright 1988-2011 The MathWorks, Inc.
%     

% Input parsing and validation
narginchk(2,6);

[method,k,Ni,traceFlag,N,NW] = parseinputs_dpss(N,NW,varargin{:});

switch method
   
case 'calc'
   if traceFlag,
      disp(getString(message('signal:dpss:ComputingTheDPSSUsingDirectAlgorithm')))
   end
   [E,V] = dpsscalc(N,NW,k);
   
case {'spline','linear'}
   if isempty(Ni)
      ind = dpssdir(NW,'NW');
      if ~isempty(ind)
         Nlist = [ind.N];
         % find closest length and use that one
         [dum,i] = min(abs(N-Nlist)); %#ok
         Ni = Nlist(i);
         if strcmp(method,'linear') && Ni < N
            if i<length(Nlist)
               Ni = Nlist(i+1);
            else
               Ni = [];
               error(message('signal:dpss:MissingTimeBandwidthProduct2', 'NW', sprintf( '%g', NW ), 'N', sprintf( '%g', N )));
            end
         end
      else
         error(message('signal:dpss:MissingTimeBandwidthProduct1', 'NW', sprintf( '%g', NW )));
      end
   end
   
   if traceFlag,
      disp([getString(message('signal:dpss:ComputingDPSSUsing')) ' ' method ' ' getString(message('signal:dpss:InterpolationFromLength')) ' ' int2str(Ni) '...'])
   end
   
   [E,V]=dpssint(N,NW,Ni,method);
   
otherwise
   error(message('signal:dpss:invalidMethod'))
   
end

%----------------------------------------------------------------------
function [E,V] = dpsscalc(N,NW,k)
%DPSSCALC Calculate slepian sequences.
%   [E,V] = dpsscalc(N,NW,k) uses tridieig() to get eigenvalues 1:k if k is 
%   a scalar, and k(1):k(2) if k is a matrix, of the sparse tridiagonal matrix.
%   It then uses inverse iteration using the exact eigenvalues on a starting 
%   vector with approximate shape, to get the eigenvectors required.  It then 
%   computes the eigenvalues V of the Toeplitz sinc matrix using a fast 
%   autocorrelation technique.

%   Authors: T. Krauss, C. Moler, E. Breitenberger
W=NW/N;
if nargin < 3
    k = min(round(2*N*W),N);
    k = max(k,1);
end
if length(k) == 1
    k = [1 k];
end

% Generate the diagonals
d=((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W);  % diagonal of B
ee=(1:N-1)'.*(N-1:-1:1)'/2;               % super diagonal of B

% Get the eigenvalues of B.
v = tridieig(d,[0; ee],N-k(2)+1,N-k(1)+1);
v = v(end:-1:1);
Lv = length(v);

%B = spdiags([[ee;0] d [0;ee]],[-1 0 1],N,N);
%I = speye(N,N);

% Compute the eigenvectors by inverse iteration with
% starting vectors of roughly the right shape.
E = zeros(N,k(2)-k(1)+1);
t = (0:N-1)'/(N-1)*pi;
lastwarn(''); msg = '';
for j = 1:Lv
   msg2= '';
   e = sin((j+k(1)-1)*t);
   e = tridisolve(ee,d-v(j),e,N);
   e = tridisolve(ee,d-v(j),e/norm(e),N);
   e = tridisolve(ee,d-v(j),e/norm(e),N);   

   [msg2,id2] = lastwarn('');
   if ~isempty(msg2)
       warning(message('signal:dpss:DPSS', j));
       lastwarn('');
       msg = msg2;
   end
   E(:,j) = e/norm(e);
end
lastwarn(msg);

d=mean(E);
for i=k(1):k(2)
   if rem(i,2)  % i is odd
     % Polarize symmetric dpss
       if d(i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
   else         % i is even
     % Polarize anti-symmetric dpss
       if E(2,i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
   end
end

% get eigenvalues of sinc matrix
%  Reference: [1] Percival & Walden, Exercise 8.1, p.390
s = [2*W; 4*W*sinc(2*W*(1:N-1)')];
q = zeros(size(E));
blksz = Lv;  % <-- set this to some small number if OUT OF MEMORY!!!
for i=1:blksz:Lv
    blkind = i:min(i+blksz-1,Lv);
    q(:,blkind) = fftfilt(E(N:-1:1,blkind),E(:,blkind));
end
V = q'*flipud(s);

% return 1 for any eigenvalues greater than 1 due to finite precision errors
V = min(V,1);
% return 0 for any eigenvalues less than 0 due to finite precision errors
V = max(V,0);

%---------------------------------------------------------------------
function [En,V] = dpssint(N, NW, M, int, E,V)
% Syntax: [En,V]=dpssint(N,NW); [En,V]=dpssint(N,NW,M,'spline');
%  Dpssint calculates discrete prolate spheroidal
%  sequences for the parameters N and NW. Note that
%  NW is normally 2, 5/2, 3, 7/2, or 4 - not i/N. The 
%  dpss are interpolated from previously calculated 
%  dpss of order M (128, 256, 512, or 1024). 256 is the 
%  default for M. The interpolation can be 'linear' 
%  or 'spline'. 'Linear' is faster, 'spline' the default.
%  Linear interpolation can only be used for M>N. 
%  Returns:
%              E: matrix of dpss (N by 2NW)
%              V: eigenvalue vector (2NW)
% 
% Errors in the interpolated dpss are very small but should be 
% checked if possible. The differences between interpolated
% values and values from dpsscalc are generally of order
% 10ee-5 or better. Spline interpolation is generally more
% accurate. Fractional errors can be very large near
% the zero-crossings but this does not seriously affect
% windowing calculations. The error can be reduced by using
% values for M which are close to N.
%
% Written by Eric Breitenberger, version date 10/3/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

if     nargin==2,
  M=256; int='spline';
elseif nargin==3,
  if ischar(M), int=M; M=256; 
  else          int='spline'; end
end

if strcmp(int, 'linear') && N > M
  error(message('signal:dpss:CannotLinearlyInterpolate', 'N', 'Ni'));
end

if nargin<=4
    [E,V] = dpssload(M,NW);
else
    if size(E,1)~=M
        error(message('signal:dpss:NiEMismatchedRowSize', 'Ni', 'E'));
    end
end

k=min(round(2*NW),N); % Return only first k values
k = max(k,1);
E=E(:,1:k);
V=V(1:k);
x=1:M;

% The scaling for the interpolation:
% This is not necessarily optimal, and 
% changing s can improve accuracy.
 
s=M/N;
midm=(M+1)/2;
midn=(N+1)/2;
delta=midm-s*midn;
xi=linspace(1-delta, M+delta, N);

% Interpolate from M values to N
% Spline interpolation is a bit better,
% but takes about twice as long.
% Errors from linear interpolation are 
% usually smaller than errors from scaling.

En=interp1(x,E,xi,['*' int]);

% Re-normalize the eigenvectors
En=En./(ones(N,1)*sqrt(sum(En.*En)));

%----------------------------------------------------------------------
function [method,k,Ni,traceFlag,N,NW] = parseinputs_dpss(N,NW,varargin)
% PARSEINPUTS Parses the inputs to the DPSS function and validates it.
%
% Inputs:   -  The inputs to this function are the same as the ones 
%              passed to DPSS.  See the help for DPSS.
%
% Outputs:
%   method    - method updated with value entered by the user
%   k         - number of most band-limited DPSSs
%   Ni        - length of DPSSs to be interpolated
%   traceFlag - a boolean flag indicating if 'trace' was resquested
%

% Here are all possible input combinations in varargin (after N,NW)...
%
% 1 Option specified
% (N,NW,k), (N,NW,method) or (N,NW,'trace')
%
% 2 Options Specified.
% (N,NW,k,method) or (N,NW,k,'trace') or  (N,NW,method,'trace')
% or (N,NW,method,Ni)
%
% 3 Options Specified.
% (N,NW,k,method,Ni), (N,NW,k,method,'trace') or (N,NW,'method',Ni,'trace')
%
% 4 Options Specified.
% (N,NW,k,method,Ni,'trace')

% Defined defaults
% If user didn't specify method or traceFlag set them to defaults.
method    = 'calc';
k         = [];
Ni        = [];
traceFlag = [];  % It gets set to 0 if user does not specified it.

% Validate input arguments N and NW.
if length(N)>1,
   error(message('signal:dpss:NeedScalarN', 'N'));
end
if length(NW)>1,
   error(message('signal:dpss:NeedScalarNW', 'NW'));
end

if ~isnumeric(N) || N < 1 || ~isfinite(N) || rem(N,1)
   error(message('signal:dpss:NeedPositiveIntegerN', 'N'));
end

if ~isnumeric(NW) || NW < 0 || ~isfinite(NW)
   error(message('signal:dpss:NeedPositiveNW', 'NW'));
end

if NW >= N/2,
   error(message('signal:dpss:AliasedNW', 'NW'));
end

N = double(N);
NW = double(NW);

% Default number of sequences to return
k = min(round(2*NW),N);
k = max(k,1);

% Validate and parse the optional input arguments
nopts = length(varargin);

for opt_indx = 1:nopts,
   arg = varargin{opt_indx};
   
   % Parse strings
   if ischar(arg),
      [method,traceFlag] = ...
         parseStrOpts(arg,method,traceFlag,nopts,opt_indx);
      
   else % Parse numerics
      
      % 1 Option Specified.
      if opt_indx == 1,
         k = arg;
         if isempty(k) || any(k ~= round(abs(k))) || any(k > N),
            error(message('signal:dpss:KOutOfRangeOfN', 'K'));
         elseif length(k) > 2 || numel(k) > 2,
            error(message('signal:dpss:BadDimensionsForK', 'K'));
         end
         k = double(k);
         
      % 2 or 3 Options Specified
      elseif any(opt_indx==[2,3]),
         Ni = arg; 
         if length(Ni) > 1 || isempty(Ni),
            error(message('signal:dpss:NeedPositiveIntegerNi', 'Ni'));
         end  
         Ni = double(Ni);
      end
      
   end 
   
end % for-loop

if isempty(traceFlag),  % If user didn't specify it set it to 0 (default).
   traceFlag = 0;
end


%----------------------------------------------------------------------
function [method,traceFlag] = ...
   parseStrOpts(inputStr,method,traceFlag,nopts,opt_indx)
% PARSESTROPTS Parses the input string options passed to dpss.
%
% Inputs:
%   inputStr  - input string entered by the user
%   method    - interpolation method
%   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace
%   nopts     - number of optional input arguments specified
%   opt_indx  - index of the input option being parsed
%
% Outputs:
%   method    - interpolation method specified by the user
%   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace

% Define output variable in case of early return.

% Define all possible string options; lower case, no punctuation:
strOpts  = {'calc','spline','linear',...
      'trace','notrace'};

i = find(strncmpi(inputStr, strOpts, length(inputStr)));
if isempty(i),
   error(message('signal:dpss:BadStringOpt', 'calc', 'spline', 'linear', 'trace'));
else
   % Determine option set to which this string applies:
   switch i
   case {1,2,3},
      method = strOpts{i};
   case 4,
      traceFlag = 1;
   case 5,         
      traceFlag = 0;
   end
end

if ~isempty(traceFlag) && nopts > opt_indx,
   error(message('signal:dpss:MisplacedTrace', 'trace')); 
end

% [EOF] dpss.m

function y = mean(x,dim,flag)
%MEAN   Average or mean value.
%   S = MEAN(X) is the mean value of the elements in X
%   if X is a vector. For matrices, S is a row
%   vector containing the mean value of each column.
%   For N-D arrays, S is the mean value of the
%   elements along the first array dimension whose size
%   does not equal 1.
%   If X is floating point, that is double or single,
%   S has the same class as X. If X is of integer
%   class, S has class double.
%
%   MEAN(X,DIM) takes the mean along the dimension DIM
%   of X.
%
%   S = MEAN(X,'double') and S = MEAN(X,DIM,'double')
%   returns S in double, even if X is single.
%
%   S = MEAN(X,'native') and S = MEAN(X,DIM,'native')
%   returns S in the same class as X.
%
%   S = MEAN(X,'default') and S = MEAN(X,DIM,'default')
%   are equivalent to S = MEAN(X) and S = MEAN(X,DIM)
%   respectively.
%
%   Example: If X = [1 2 3; 3 3 6; 4 6 8; 4 7 7];
%
%   then mean(X,1) is [3.0000 4.5000 6.0000] and
%   mean(X,2) is [2.0000 4.0000 6.0000 6.0000].'
%
%   Class support for input X:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32,
%               int32, uint64, int64
%
%   See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.

%   Copyright 1984-2013 The MathWorks, Inc.

if nargin==2 && ischar(dim)
    flag = dim;
elseif nargin < 3
    flag = 'default';
end

if nargin == 1 || (nargin == 2 && ischar(dim))
    % preserve backward compatibility with 0x0 empty
    if isequal(x,[])
        y = sum(x,flag)/0;
        return
    end
    
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

if ~isobject(x) && isinteger(x) 
    isnative = strcmp(flag,'native');
    if intmin(class(x)) == 0  % unsigned integers
        y = sum(x,dim,flag);
        if (isnative && all(y(:) < intmax(class(x)))) || ...
                (~isnative && all(y(:) <= flintmax))
            % no precision lost, can use the sum result
            y = y/size(x,dim);
        else  % throw away and recompute
            y = intmean(x,dim,flag);
        end
    else  % signed integers
        ypos = sum(max(x,0),dim,flag);
        yneg = sum(min(x,0),dim,flag);
        if (isnative && all(ypos(:) < intmax(class(x))) && ...
                all(yneg(:) > intmin(class(x)))) || ...
                (~isnative && all(ypos(:) <= flintmax) && ...
                all(yneg(:) >= -flintmax))
            % no precision lost, can use the sum result
            y = (ypos+yneg)/size(x,dim);
        else  % throw away and recompute
            y = intmean(x,dim,flag);    
        end
    end
else
    y = sum(x,dim)/size(x,dim);
end

function y = intmean(x, dim, flag)
% compute the mean of integer vector

if strcmp(flag, 'native')
    doubleoutput = false;
else 
    assert(strcmp(flag,'double') || strcmp(flag,'default'));
    doubleoutput = true;
end

shift = [dim:ndims(x),1:dim-1];
x = permute(x,shift);

xclass = class(x);
if doubleoutput
    outclass = 'double';
else
    outclass = xclass;
end

if intmin(xclass) == 0
    accumclass = 'uint64';
else
    accumclass = 'int64';
end
xsiz = size(x);
xlen = cast(xsiz(1),accumclass);

y = zeros([1 xsiz(2:end)],outclass);
ncolumns = prod(xsiz(2:end));
int64input = isa(x,'uint64') || isa(x,'int64');

for iter = 1:ncolumns
    xcol = cast(x(:,iter),accumclass);
    if int64input
        xr = rem(xcol,xlen);
        ya = sum((xcol-xr)./xlen,1,'native');
        xcol = xr;
    else
        ya = zeros(accumclass);
    end
    xcs = cumsum(xcol);
    ind = find(xcs == intmax(accumclass) | (xcs == intmin(accumclass) & (xcs < 0)) , 1);
    
    while (~isempty(ind))
        remain = rem(xcs(ind-1),xlen);
        ya = ya + (xcs(ind-1) - remain)./xlen;
        xcol = [remain; xcol(ind:end)];
        xcs = cumsum(xcol);
        ind = find(xcs == intmax(accumclass) | (xcs == intmin(accumclass) & (xcs < 0)), 1);
    end
    
    if doubleoutput
        remain = rem(xcs(end),xlen);
        ya = ya + (xcs(end) - remain)./xlen;
        % The latter two conversions to double never lose precision as
        % values are less than FLINTMAX. The first conversion may lose
        % precision.
        y(iter) = double(ya) + double(remain)./double(xlen);
    else
        y(iter) = cast(ya + xcs(end) ./ xlen, outclass);
    end
end
y = ipermute(y,shift);

%EOF

function y = fftfilt(b,x,nfft)
%FFTFILT Overlap-add method for FIR filtering using FFT.
%   Y = FFTFILT(B,X) filters X, with the FIR filter specified by the vector
%   of coefficients B, using the overlap/add method, and internal
%   parameters (FFT size and block length) that guarantee efficient
%   execution.
%   
%   Y = FFTFILT(B,X,N) allows you to have some control over the internal
%   parameters, by using an FFT of at least N points.
%
%   Y = FFTFILT(D,X,...) filters X with the FIR digital filter D. You
%   design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   If X is a matrix, FFTFILT filters its columns.  If B is a matrix,
%   FFTFILT applies the filter in each column of B to the signal vector X.
%   If B and X are both matrices with the same number of columns, then the
%   i-th column of B is used to filter the i-th column of X.
%
%   It is advantageous to use FFTFILT instead of FILTER when the signal is
%   relatively large.  FILTER performs N multiplications for each sample in
%   X where N is the filter length.  FFTFILT performs 2 FFT operations at
%   the cost of L*log2(L)/2 where L is the block length.  It then performs
%   L pointwise multiplications for a total cost of L*(1+log2(L))
%   multiplications.  The cost ratio is therefore L*(1+log2(L))/(N*L) =>
%   (1+log2(L))/N which is approximately log2(L)/N.  Therefore FFTFILT
%   becomes advantageous when log2(L) is less than N.
%
%   % Example 1:
%   %   Construct a Signal and filter it with a 10 point averaging filter
%   %   using fftfilt. 
%
%   fs = 100;                               % Sampling frequency
%   t = 0:1/fs:1;                           % Time vector
%   x = sin(2*pi*t*3)+.25*sin(2*pi*t*40);   % Input Signal
%   b = ones(1,10)/10;  % 10 point averaging filter
%   y = fftfilt(b,x);   % FIR filtering using overlap-add method
%   plot(t,x,t,y,'--');
%   legend('Original Signal','Filtered Signal')
%
%   % Example 2:
%   %   Use the designfilt function to design a lowpass FIR digital filter 
%   %   with order 350 and cutoff frequency of 150 Hz. The sample rate is 
%   %   1.5 KHz. Filter a long vector of data using the overlap-add method 
%   %   to increase speed.
%
%   D = designfilt('lowpassfir', 'FilterOrder', 350, ...
%    'CutoffFrequency', 150, 'SampleRate', 1500);
%
%   data = randn(10e6,1);  
%   y = fftfilt(D,data);
%
%   See also FILTER, FILTFILT.

%   --- Algorithmic details ---
%   The overlap/add algorithm convolves B with blocks of X, and adds
%   the overlapping output blocks.  It uses the FFT to compute the
%   convolution.
% 
%   Particularly for long FIR filters and long signals, this algorithm is 
%   MUCH faster than the equivalent numeric function FILTER(B,1,X).
%
%   Y = FFTFILT(B,X) -- If you leave N unspecified:   (RECOMMENDED)
%       Usually, length(X) > length(B).  Here, FFTFILT chooses an FFT 
%       length (N) and block length (L) which minimize the number of 
%       flops required for a length-N FFT times the number of blocks
%       ceil(length(X)/L).  
%       If length(X) <= length(B), FFTFILT uses a single FFT of length
%       nfft = 2^nextpow2(length(B)+length(X)-1), essentially computing 
%       ifft(fft(B,nfft).*fft(X,nfft)).
%
%   Y = FFTFILT(B,X,N) -- If you specify N:
%       In this case, N must be at least length(B); if it isn't, FFTFILT 
%       sets N to length(B).  Then, FFTFILT uses an FFT of length 
%       nfft = 2^nextpow2(N), and block length L = nfft - length(B) + 1. 
%       CAUTION: this can be VERY inefficient, if L ends up being small.

%   Author(s): L. Shure, 7-27-88
%              L. Shure, 4-25-90, revised
%              T. Krauss, 1-14-94, revised
%   Copyright 1988-2013 The MathWorks, Inc.

%   Reference:
%      A.V. Oppenheim and R.W. Schafer, Digital Signal 
%      Processing, Prentice-Hall, 1975.

narginchk(2,3);

m = size(x, 1);
if m == 1
    x = x(:);    % turn row into a column
end

nx = size(x,1);

if min(size(b))>1
   if (size(b,2)~=size(x,2))&&(size(x,2)>1)
      error(message('signal:fftfilt:InvalidDimensions'))
   end
else
   b = b(:);   % make input a column
end
nb = size(b,1);

if nargin < 3
% figure out which nfft and L to use
    if nb >= nx || nb > 2^20    % take a single FFT in this case
        nfft = 2^nextpow2(nb+nx-1);
        L = nx;
    else
        fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
       149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];
        n = 2.^(1:20); 
        validset = find(n>(nb-1));   % must have nfft > (nb-1)
        n = n(validset); 
        fftflops = fftflops(validset);
        % minimize (number of blocks) * (number of flops per fft)
        L = n - (nb - 1);
        [dum,ind] = min( ceil(nx./L) .* fftflops ); %#ok
        nfft = n(ind);
        L = L(ind);
    end

else  % nfft is given
  % Cast to enforce precision rules
    nfft = signal.internal.sigcasttofloat(nfft,'double','fftfilt','N','allownumeric');
    if nfft < nb
        nfft = nb;
    end
    nfft = 2.^(ceil(log(nfft)/log(2))); % force this to a power of 2 for speed
    L = nfft - nb + 1;
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(b,x,nfft);
catch ME
    throwAsCaller(ME);
end

B = fft(b,nfft);
if length(b)==1,
     B = B(:);  % make sure fft of B is a column (might be a row if b is scalar)
end
if size(b,2)==1
    B = B(:,ones(1,size(x,2)));  % replicate the column B 
end
if size(x,2)==1
    x = x(:,ones(1,size(b,2)));  % replicate the column x 
end
y = zeros(size(x));

istart = 1;
while istart <= nx
    iend = min(istart+L-1,nx);
    if (iend - istart) == 0
        X = x(istart(ones(nfft,1)),:);  % need to fft a scalar
    else
        X = fft(x(istart:iend,:),nfft);
    end
    Y = ifft(X.*B);
    yend = min(nx,istart+nfft-1);
    y(istart:yend,:) = y(istart:yend,:) + Y(1:(yend-istart+1),:);
    istart = istart + L;
end

if ~(any(imag(b(:))) || any(imag(x(:))))
	y = real(y);
end

if (m == 1)&&(size(y,2) == 1)
    y = y(:).';    % turn column back into a row
end

%EOF
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

%EOF

function y = flipud(x)
%FLIPUD Flip matrix in up/down direction.
%   FLIPUD(X) returns X with columns preserved and rows flipped
%   in the up/down direction.  For example,
%   
%   X = 1 4      becomes  3 6
%       2 5               2 5
%       3 6               1 4
%
%   Class support for input X:
%      float: double, single
%
%   See also FLIPLR, ROT90, FLIPDIM.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.9.4.5 $  $Date: 2010/08/23 23:08:00 $

if ~ismatrix(x)
  error(message('MATLAB:flipud:SizeX')); 
end
y = x(end:-1:1,:);

%EOF

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

%EOF

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
% [EOF]

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

function y=sinc(x)

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   

function chkinputdatatype(varargin)
%CHKINPUTDATATYPE Check that all inputs are double

%   Copyright 2009-2013 The MathWorks, Inc.


for n = 1:nargin
    if ~isa(varargin{n},'double')
        error(message('signal:chkinputdatatype:NotSupported'));
    end
end



% [EOF]

function x = tridieig(c,b,m1,m2,eps1);
%TRIDIEIG  Find a few eigenvalues of a tridiagonal matrix.
%   LAMBDA = TRIDIEIG(D,E,M1,M2).  D and E, two vectors of length N,
%   define a symmetric, tridiagonal matrix:
%      A = diag(E(2:N),-1) + diag(D,0) + diag(E(2:N),1)
%   E(1) is ignored.
%   TRIDIEIG(D,E,M1,M2) computes the eigenvalues of A with indices
%      M1 <= K <= M2.
%   TRIDIEIG(D,E,M1,M2,TOL) uses TOL as a tolerance.
%
%   See also TRIDISOLVE.

%   Author: C. Moler
%   Copyright 1988-2002 The MathWorks, Inc.

%   Mex-file version will be called if on the path.

%mbrealvector(c);
%mbrealvector(b);
%mbintscalar(m1);
%mbintscalar(m2);
%mbrealscalar(eps1)
if nargin < 5, eps1 = 0; end
n = length(c);
b(1) = 0;
beta = b.*b;
xmin = min(c(n) - abs(b(n)),min(c(1:n-1) - abs(b(1:n-1)) - abs(b(2:n))));
xmax = max(c(n) + abs(b(n)),max(c(1:n-1) + abs(b(1:n-1)) + abs(b(2:n))));
eps2 = eps*max(xmax,-xmin);
if eps1 <= 0, eps1 = eps2; end
eps2 = 0.5*eps1 + 7*eps2;

x0 = xmax;
x = zeros(n,1);
wu = zeros(n,1);
x(m1:m2) = xmax(ones(m2-m1+1,1));
wu(m1:m2) = xmin(ones(m2-m1+1,1));
z = 0;
for k = m2:-1:m1
   xu = xmin;
   for i = k:-1:m1
      if xu < wu(i)
         xu = wu(i);
         break
      end
   end
   if x0 > x(k), x0 = x(k); end
   while 1
      x1 = (xu + x0)/2;
      if x0 - xu <= 2*eps*(abs(xu)+abs(x0)) + eps1
         break
      end
      z = z + 1;
      a = 0;
      q = 1;
      for i = 1:n
         if q ~= 0
            s = beta(i)/q;
         else
            s = abs(b(i))/eps;
         end
         q = c(i) - x1 - s;
         a = a + (q < 0);
      end
      if a < k
         if a < m1
            xu = x1;
            wu(m1) = x1;
         else
            xu = x1;
            wu(a+1) = x1;
            if x(a) > x1, x(a) = x1; end
         end
      else
         x0 = x1;
      end
   end
   x(k) = (x0 + xu)/2;
end
x = x(m1:m2)';

%TRIDISOLVE  Solve A*x = b where A is a square, symmetric tridiagonal matrix.
%   X = TRIDISOLVE(E,D,B,N) is the solution to the system
%     A*X = B
%   where A is an N-by-N symmetric, real, tridiagonal matrix given by
%     A = diag(E,-1) + diag(D,0) + diag(E,1)
%   Algorithm from Golub and Van Loan, "Matrix Computations", 2nd Edition, 
%   p.156.
%   Assumes A is non-singular.  If A is singular, a warning is issued and
%   results may be inaccurate.
%
%   Called by DPSS.
%
%   See also TRIDIEIG.

%   MEX-file implementation: T. Krauss, D. Orofino, 4/22/97
%   Copyright 1988-2002 The MathWorks, Inc.

%   MEX-file.


