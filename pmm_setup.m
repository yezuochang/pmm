function setup_pmm()
setup_path;
setup_cvx;
%setup_pso;
preinstall_methods;

function preinstall_methods()
pmm_install('VF', 'Vector Fitting', 0, 0);
pmm_install('SDP', 'SDP Method', 0, 1);
pmm_install('LC', 'Local Compensation', 0, 1);
pmm_install('EPM', 'Eigenvalue Perturbation', 0, 1);
pmm_install('FRP', 'FRP Method', 0, 1);
pmm_install('DAO', 'Domain Alternated Optimization', 1, 0);
pmm_install('EPM2', 'Eigenvalue Perturbation improved', 0, 1);
pmm_install('EPM3', 'Eigenvalue Perturbation improved', 0, 1);
pmm_install('EPM4', 'Eigenvalue Perturbation improved', 0, 1);

function setup_path()
[root,name,ext]=fileparts(mfilename('fullpath'));
addpath(root);
addpath(sprintf('%s/utils',root));
addpath(sprintf('%s/methods',root));
addpath(sprintf('%s/methods/dao',root));
addpath(sprintf('%s/methods/epm',root));
addpath(sprintf('%s/others/mfit',root));
addpath(sprintf('%s/others/cvx',root));

function setup_pso()
[root,name,ext]=fileparts(mfilename('fullpath'));
dir=pwd;
cd(sprintf('%s/others/PSOt',root));
pso_setup;
cd(dir);
cvx_loaded=1;


function setup_cvx()
   
[root,name,ext]=fileparts(mfilename('fullpath'));
dir=pwd;
cd(sprintf('%s/others/cvx',root));
cvx_setup;
cd(dir);
cvx_loaded=1;
