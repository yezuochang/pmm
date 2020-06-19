function [F,H]=ldstone(filename,opts);
% ldstone -- load touchstone format data
% Input 
%     filename: filename of sNp data file.
%     opts    : options
% Output
%     F       : A vector of length Ns containing the frequency samples.
%     H       : An array of size m-by-m-by-Ns containing the data
% Options (opts)    
%     opts.parametertype = 
if nargin < 2
    opts=pmm_default;
end
% data=read(rfdata.data,filename);
% F=data.Freq;
% S=data.S_Parameters;
% [F,S]=SXPParse(filename);



sobj = sparameters(filename);
F = sobj.Frequencies;
S = sobj.Parameters;

% remove samples with freq>fmax
fmax = optget(opts,'fmax',1e100);
ix=find(F>fmax);
F=F(1:ix);
S=S(:,:,1:ix);

if ~opts.enforceDC
    avg=sum(F)/length(F);
    ix=find(F<avg*1e-10);    
    if ~isempty(ix)
        fprintf('Warning: The following frequencies are removed as they are too small comparing to other frequencies. \n');
        fprintf('%e\n', F(ix));
        px=setdiff(1:length(F),ix);
        F=F(px);
        S=S(:,:,px);
    end
end

if optget(opts,'parametertype','Y') == 'Y'
    H=s2y(S,50);
else
    H=S;
end

ns=size(H,3);
for c=1:ns
    H1(:,:,c)=(H(:,:,c)+H(:,:,c).')/2;
end
H=H1;
