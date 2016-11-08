% Example of continuation power flow
close all;clc;clear;
dirCPF = ones(3,1); % Load increase at the 3 load buses
results_cpf = ch_runCPF('case9static','loads579',0,dirCPF);
processResCPF('case9static','loads579');
plotAfterCPF('case9static','loads579');