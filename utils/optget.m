function value=optget(opts,name,def)
if ~isa(opts,'struct')
    value = def;
    return;
end

if isfield(opts,name)
    value=opts.(name);
else
    value=def;
end