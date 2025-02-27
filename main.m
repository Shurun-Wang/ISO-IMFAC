%___________________________________________________________________%
%  Snake Optimizer (SO) source codes version 1.0                    %
%                                                                   %
%  Developed in MATLAB R2021b                                       %
%                                                                   %
%  Author and programmer:  Fatma Hashim & Abdelazim G. Hussien      %
%                                                                   %
%         e-Mail: fatma_hashim@h-eng.helwan.edu.eg                  %
%                 abdelazim.hussien@liu.se                          %
%                 aga08@fayoum.edu.eg                               %
%                                                                   %
%                                                                   %
%   Main paper: Fatma Hashim & Abdelazim G. Hussien                 %
%               Knowledge-based Systems                             %
%               in press,                                           %
%               DOI: 10.1016/j.knosys.2022.108320                   %
%                                                                   %
%___________________________________________________________________%
close all
clear all
clc
% fix the random seed (random select)

rng(10)
fitfun = @IMFAC;
dim=8;
Max_iteration=100;
SearchAgents_no = 40;
ub = [1 1 20 20 20 100 100 10];
lb = [1e-7 1e-7 1e-7 1e-7 1e-7 0 0 0];
tlt='IMFAC';
tic
[Xfood, Xvalue,CNVG] = ISO(SearchAgents_no,Max_iteration,fitfun, dim,lb,ub);
toc
Xvalue
Xfood

% hold on
% plot(CNVG,'Color', 'r')
% xlim([1 100]);
