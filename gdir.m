function [newdir] = gdir
% gdir.m
%
% simple code to change to a new working directory
% uses matlab GUI

newdir = uigetdir;
cd(newdir)

