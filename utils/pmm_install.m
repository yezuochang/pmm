function pmm_install(funcname, description, passivity_in, passivity_out)
global pmm_methods;
pmm_methods.(funcname).func=eval(sprintf('@%s',funcname));
pmm_methods.(funcname).description=description;
pmm_methods.(funcname).passive_in=passivity_in;
pmm_methods.(funcname).passive_out=passivity_out;
