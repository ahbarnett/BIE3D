function [T,Tt] = pulsegaussian(a0,t0,s0)
% PULSEGAUSSIAN  function handles for pulse func of time
%
% [T,Tt] = pulsegaussian(a0,t0,s0) returns function handle for T(t) and T'(t)
%  of a Gaussian pulse of width s0, centered at time t0, with amplitude a0.
T = @(t) a0*exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
