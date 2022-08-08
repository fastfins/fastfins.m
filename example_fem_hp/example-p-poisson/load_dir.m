
root = pwd();
%
fin_path = [root '/../..'];
warning('off')
rmpath(genpath(fin_path));
warning('on')
addpath(genpath([fin_path '/fastfins/']));
addpath(genpath([fin_path '/fem_hp']));
addpath(genpath([fin_path '/rom_tools']));
%
addpath(genpath([root '/model']));
addpath(genpath([root '/rom']));
%
fastfins_check_solvers()

