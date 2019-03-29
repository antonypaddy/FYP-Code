function [out] = resize(A,si)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

out = [A; zeros(si - length(A),1)];

end

