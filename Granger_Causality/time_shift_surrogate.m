function xs = time_shift_surrogate(x, minShift)
% TIME_SHIFT_SURROGATE  Circular time shift surrogate
%
% Inputs:
%   x        : T x 1 vector
%   minShift: minimum shift (e.g. 20)
%
% Output:
%   xs : surrogate signal

    T = length(x);

    l = randi([minShift, T-minShift]);

    xs = [x(l+1:end); x(1:l)];
end
