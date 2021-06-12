function [N, bin_centers] = histbins(data, bin_edges)
% HISTBINS - Generate counts and bin centers for custom histogram
% 
%  [N, BIN_CENTERS] = HISTBINS(DATA, BIN_EDGES)
%
%   Inputs: DATA - 1-dimensional data observations
%           BIN_EDGES - The edges of the bins
%   Outputs: N - A vector of counts for each bin
%            BIN_CENTERS - The bin centers

N = histc(data, bin_edges);
bin_centers = 0.5*(bin_edges(2:end) + bin_edges(1:end-1));
N = N(1:end-1);