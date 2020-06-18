function [G,W,F,H,info]=pmm(inputfile,q,opts)
% pmm -- Passive Macro Modeling
% Please make sure to run pmm_setup.m first.
% 
% Usage: 
% [G, W, F, H] = pmm(infile, q, opts)
% See demo/pmm_demo.m
% 
% Copyright 2012, Zuochang Ye, zuochang@tsinghua.edu.cn
%
if nargin < 3
    opts = [];
	opts=pmm_default;
end
opts.q=q;

[inputfile,outputfile,subcktname]=ParseInput(inputfile);

%% Step 0: load data
[F,H]=ldstone(inputfile,opts);

% [F,H]=freqinterp(F,H,opts);

opts=select_method(size(H,1),q,opts);
check_methods(opts);
opts
[F1,H1,scale]=datascale(F,H,opts);
G=[]; W=[];
if isfield(opts,'Func')
    fprintf('\n');
    info=[];
    for k=1:length(opts.Func)
        [G,W,info_k]=func_call(opts.Func{k},G,W,F1,H1,opts);
        info{k}=info_k;
        if ~info_k.success
           break;
        end
    end
end
% if (isempty(W) || optget(opts,'sdp',0))
%     W=spf(G);
% end

%% Step 4: scale back
G=ss_scale(G,scale);
if exist('W','var') && ~isempty(W)
    W=ss_scale(W,scale);

    ss_export(G,W,subcktname,outputfile,opts);
    ss_verify(G,W,F,H,opts);
end
print_info(info);

function print_info(info)
fprintf('  FuncName\t     Time     \t     Error     \tPassivity\n');
fprintf('------------------------------------------------------\n');
for c=1:length(info)
    fprintf('%10s\t%e\t%e\t%s\n',...
        info{c}.func,info{c}.time,info{c}.error,info{c}.passivity);
end

function [G,W,info]=func_call(funcname,G,W,F1,H1,opts)
global pmm_methods;

func        = pmm_methods.(funcname).func;
passive_in  = pmm_methods.(funcname).passive_in;
passive_out = pmm_methods.(funcname).passive_out;

t0=cputime;
[G,W]=func(G,W,F1,H1,opts);
t=cputime-t0;

r2=passivity_violation(G);
if passive_out && ~isempty(r2)
    fprintf('Warning: %s claims but fails to ensure passivity of the system.\n', funcname);
    fprintf('Check the following frequencies for diagnosis.\n');
    fprintf('%e\n',r2/2/pi);
    info.success = 0;
else
    info.success=1;
end

wgt_scheme=optget(opts,'wgt_scheme',1);
info.func=func2str(func);
info.time=t;
info.error=ss_error(G,F1,H1,wgt_scheme);
if isempty(r2)
    info.passivity='passive';
else
    info.passivity='non-passive';
end

%fprintf('%10s\t%e\t%e\t%d\n',...
%    info{c}.func,info{c}.time,info{c}.error,info{c}.passivity);

%% ParseInput()
function [inputfile,outputfile,subcktname]=ParseInput(inputfile)
[filepath,filename,ext]=fileparts(inputfile);
if isempty(filepath)
	filepath='.';
end
outputfile=sprintf('%s/%s_ncss.scs',filepath,filename);
subcktname=sprintf('%s_ncss',filename);



function opts=select_method(m, q, opts)

method=optget(opts,'method','auto');
if strcmp(method,'user_defined')
    return;
end
opts.Func={};
opts.Func{1} = 'VF';
if strcmp(method, 'auto')
    n=m*q;
    if n*m < 100
        opts.Func{2} = 'SDP';
    else        
        opts.Func{2} = 'LC';
        opts.Func{3} = 'DAO';        
    end
    return;
end

if strcmp(method,'sdp')
    opts.Func{2} = 'SDP';
    return;
end

if strcmp(method,'vf_only')
    return;
end

if strcmp(method, 'lc_only')
    opts.Func{2} = 'LC';    
    return;
end

if strcmp(method, 'epm_only')
    opts.Func{2} = 'EPM';    
    return;
end

if strcmp(method, 'frp_only')
    opts.Func{2} = 'FRP';    
    return;
end

if strcmp(method, 'lc_dao')
    opts.Func{2} = 'LC';
    opts.Func{3} = 'DAO';
    %opts.Func{4} = 'EPM';
    return;
end

if strcmp(method, 'epm_dao')
    opts.Func{2} = 'EPM';    
    opts.Func{3} = 'DAO';
    %opts.Func{4} = 'EPM';
    return;
end

if strcmp(method, 'frp_dao')
    opts.Func{2} = 'FRP'; 
    opts.Func{3} = 'DAO';
    %opts.Func{4} = 'EPM';
    return;
end

if strcmp(method, 'epm2_only')
    opts.Func{2} = 'EPM2';    
    return;
end

if strcmp(method, 'epm3_only')
    opts.Func{2} = 'EPM3';    
    return;
end

if strcmp(method, 'epm4_only')
    opts.Func{2} = 'EPM4';    
    return;
end

function check_methods(opts)
global pmm_methods;

for c=1:length(opts.Func)
    if ~isfield(pmm_methods,opts.Func{c})
        error(sprintf('Method %s not installed.', opts.Func{c}));
    end
end
