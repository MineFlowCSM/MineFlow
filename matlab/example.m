%Run this first! (but only once) to create the .mex file in your project
%mex COMPFLAGS='$COMPFLAGS /std:c++17 /I "..\"' "mineflowmex.cpp"

% the data files are x fastest, then y, then z. So must transpose 2d data
% files for matlab's column major ordering
ebv = int64(transpose(reshape(readmatrix("../data/sim2d76.dat"), [75, 40])));
upit = mineflowmex(ebv, 45);
value = sum(times(ebv, int64(upit)),"all");
% value should be 295932

% and must permute the 3d ones.
ebv = int64(permute(reshape(readmatrix("../data/bauxitemed.dat"), [120,120,26]), [3, 2, 1]));
upit = mineflowmex(ebv, 45);
value = sum(times(ebv, int64(upit)),"all");
% value should be 28288679

