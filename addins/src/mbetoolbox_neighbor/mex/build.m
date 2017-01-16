if (isunix)
	mex -Dfalse=0 -Dtrue=1 neighbor_mex.c phylip.c dist.c
else
	mex -DWIN32 -Dfalse=0 -Dtrue=1 neighbor_mex.c phylip.c dist.c
    %mex -DWIN64 -Dfalse=0 -Dtrue=1 neighbor_mex.c phylip.c dist.c
end

