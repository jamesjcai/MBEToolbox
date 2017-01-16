function [model] = model_nt(name,varargin)
%MODEL_NT - Returns a structure of nucleotide substitution model
%
% Syntax: [model] = model_nt('name','rmatrix',rmatrix,'freq',freq)
%
% Inputs:
%    name      - {'jc'|'k2p'|'f81'|'f84'|'hky'|'gtr'}
%    rmatrix   - Rate matrix parameters, 1x5 or 1x6 
%                e.g. [1.0 1.33333 1.0 1.0 1.333333 1.0]
%    freq      - Equilibrium frequency parameters, 1x3 or 1x4 
%                e.g. [.1 .2 .3 .4]
%
% Outputs:
%    model  - model strucure
%
% See also: 

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

if (nargin==0)
	name='jc';
end

modelNames = {'jc','k2p','f81','f84','hky','gtr','rev'};
modelName = strmatch(lower(name),modelNames); %#ok

  if isempty(modelName) 
      error('MBEToolbox:NotValidModel',...
            'Not a valid model of DNA substitution')
  elseif numel(modelName)>1
      error('MBEToolbox:AmbiguousModel',...
            'Ambiguous model of DNA substitution')
  end



[i_rmatrix,i_freq] = optionalArg(varargin);


%R=i_equalrmatrix;
%freq=i_equalbasefreq;
%persistent model
%if isempty model
%	model.R = i_equalrmatrix;
%	model.freq = i_equalbasefreq;
%end

%switch nargin
%	case 0 % Return the entire structure
%		name='jc';
%	case 1 % Return the parameter requested
%	case 2 % Set the parameter requested
%	S.(varargin{1}) = varargin{2};
%	% For MATLAB 6.1 and earlier, use
%	% S = setfield(S,varargin{1},varargin{2});
%end

%if (nargin<3), freq=[.25 .25 .25 .25]; end
%if (nargin<2), rmatrix=[1 1 1 1 1 1]; end
%if (nargin<1), name='jc'; end


switch (lower(name))
    case ('jc')
	R=i_equalrmatrix;
	freq=i_equalbasefreq;
    case ('k2p')
	R=r2R(i_rmatrix);
	freq=i_equalbasefreq;
    case ('f81')
	R=i_equalrmatrix;
	freq=i_confirmfreq(i_freq);
    case ('f84')
	R=r2R(i_rmatrix);
	freq=i_confirmfreq(i_freq);
%The F84 model (Felsenstein and Churchill, 1996), as implemented in DNAML 
%in the PHYLIP package (Felsenstein, 1993), is very similar to HKY but differs 
%slightly in how it treats transitions and transversions. This model requires 
%the same parameters as HKY. 
    case ('hky')
%The Hasegawa, Kishino and Yano (HKY) model (Hasegawa et al. 1985) allows for 
%a different rate of transitions and transversions as well as unequal frequencies 
%of the four nucleotides (base frequencies). The parameters requires by this model 
%are transition to transversion ratio (TS/TV) and the base frequencies. There are 
%a number of simpler models that are specific cases of the HKY modesl. If the base 
%frequencies are set equal then the model becomes equivalent to the Kimura 2-parameter 
%(K2P) model (Kimura, 1980). If the TS/TV is set to 0.5 as well, then it becomes 
%equivalent to the Jukes-Cantor (JC69) model (Jukes and Cantor, 1969). If the TS/TV is 
%et to 0.5 and the base frequencies are not equal then the model is equivalent to the 
%F81 model (Felsenstein, 1981). 
	R=r2R(i_rmatrix);
	freq=i_confirmfreq(i_freq);
    case {'rev','gtr'}
	R=r2R(i_rmatrix);	
	freq=i_confirmfreq(i_freq);
end

model.name=name;
model.R=R;
model.freq=freq;




function [R] = i_equalrmatrix
	R=zeros(4,4);
	for i=1:4, R(i,i)=0; for j=1:4, if i~=j, R(i,j)=1./3; end, end, end

function [freq] = i_equalbasefreq
	freq=zeros(1,4);
	for i=1:4, freq(i)=1./4.; end

function [freq] = i_confirmfreq(freq)
	errortag=0;
	[n,m] = size(freq);
	if ~(n==1&&m==4), 
		errortag=1;
	else
	 if (abs(sum(freq)-1)>eps)
		errortag=1;	 
	 end	 	
	end
	if (errortag==1)
		warning('freq supplied conatins error! MBEToolbox uses equal base freq instead.');
		freq=i_equalbasefreq;
	end	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmatrix,freq,kappa,R] = optionalArg(varargin)
%         
% Checks optional input parameter list for MODEL and
% returns user selected optional parameters or default values
% that will be evaluated. Any number of optional 
% key/parameter pairs can be supplied in any order, from 0 to the maximum 
% number supported (currently 6).
%
% USAGE:  MODEL_NT(name,'key1','strValue','key2',[n1,n2],'key3','str'...
%
% Here's a list of the optional parameter keys, description of parameter, 
% default value in parenthesis, and data type of parameter.  Remember, they
% can be supplied in any order, from 0 to all 4 pairs.
%
% 'freq'-4 element vector: [A, C, G, T] (number). 
% 'rmatrix'-Coast outline color.  Any valid matlab color (string).  Default 
%             is 'g' (green).
% 'R'    -
% 'kappa'    -
%

% Unpack the varargin cell array: no matter how may optional parameters are
% supplied, varargin is sent into this function as a 1 X 1 cell array 
% containing a 1 X N cell array of the optional parameters.  Use comma separated
% syntax format, varargin{:}, to extract list and reassign.

% COMMENT OUT varargin{:} when testing function with tVar.m
varargin = varargin{:};

% Make varargin like an associative array, row 1 key, and corresponding
% row 2 the value.  NOTE: varargin is a cell array.
col = size(varargin,2)/2;

% Handle case where user chooses no optional arguments.  Create a 1x2 empty cell
% array.  This allows function to run to completion and select all default 
% values.
if col < 1
   varargin = cell(1,2);
   col = 1;
end

varargin = reshape(varargin,2,col);

% Make all key values lower case in case user capitalized any.
key = char(varargin{1,:});

% strjust all keys to left in case user inputs for ex. ' times' instead 
% of 'times'.  User inputs 'times ' is OK.  strmatch doesn't like ' times'.
key = strjust(lower(key),'left');


% Search for optional arguments, if not found use default values.

% Determine the axis limits for the plot.
index = strmatch('rmatrix',key);
if isempty(index)
   rmatrix = [1 1 1 1 1 1];
else
   rmatrix = varargin{2,index};
end

index = strmatch('freq',key);
if isempty(index)
   freq = i_equalbasefreq;
else
   freq = varargin{2,index};
end

%index = strmatch('kappa',key);
%if isempty(index)
%   kappa = 4;
%else
%   kappa = varargin{2,index};
%end

%index = strmatch('R',key);
%if isempty(index)
%   R = r2R([1 1 1 1 1 1]);
%else
%   R = varargin{2,index};
%end

return


% The Hasegawa, Kishino and Yano (HKY) model (Hasegawa et al., 1985) allows for a different rate of transitions and transversions as well as unequal frequencies of the four nucleotides (base frequencies). The parameters requires by this model are transition to transversion ratio (TS/TV) and the base frequencies. There are a number of simpler models that are specific cases of the HKY model. If the base frequencies are set equal then the model becomes equivalent to the Kimura 2-parameter (K2P) model (Kimura, 1980). If the TS/TV is set to 0.5 as well, then it becomes equivalent to the Jukes-Cantor (JC69) model (Jukes and Cantor, 1969). If the TS/TV is set to 0.5 and the base frequencies are not equal then the model is equivalent to the F81 model (Feålsenstein, 1981). 
% The F84 model (Felsenstein and Churchill, 1996), as implemented in DNAML in the PHYLIP package (Felsenstein, 1993), is very similar to HKY but differs slightly in how it treats transitions and transversions. This model requires the same parameters as HKY. 
% Finally, the general reversible process (REV) model (e.g. Yang, 1994) allows 6 different rate parameters and is the most general model that is still consistent with the requirement of being reversible. The 6 parameters are the relative rates for each type of substitution (i.e. A to C, A to G, A to T, C to G, C to T and G to T). As this is a time-reversible process, the rate parameter of one type of substitution (e.g., A to T) is assumed to be the same as the reverse (e.g., T to A). 
