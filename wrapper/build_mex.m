% build mex wrapper
clc;clear;close all;
% cd to the path of where EPTAM_mapping.cp is at.
cd([pwd,'/wrapper']);

% CHANGE THIS PATH TO THE INCLUDE FOLDER OF MEXPLUS ON YOU MACHINE.
% https://github.com/kyamagu/mexplus
mexplus_path = '-I/home/zhouyi/workspace/project/mexplus/include/';

% path of eigen3
eigen_path = '-I/usr/include/eigen3/';

% compile (You may need to switch the version of gcc to 4.9.x for compatibility.)
mex('-v', 'GCC=/usr/bin/gcc-4.9',mexplus_path,eigen_path,'EPTAM_mapping.cpp', '-output', 'EPTAM_mapping_mex');
