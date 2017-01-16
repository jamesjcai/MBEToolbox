function [snpOffset,subjectID,allele1,allele2] = prettybaseread(filename)
%PRETTYBASEREAD read a Prettybase formatted file.
%
%   S = PRETTYBASEREAD(FILENAME) reads a PRETTYBASE format file FILENAME,
%   returning the data in the file as a structure. 
%
%   [OFFSET,SUBJECT,ALLELE1,ALLELE2] = PRETTYBASEREAD(FILENAME) reads the
%   file into separate variables. 
%
%   Prettybase is the format used by the SeattleSNPs Programs for Genomic
%   Applications (PGA) http://pga.gs.washington.edu .
%
%   Example:
%
%       % Read the sequence for the human p53 tumor gene.
%       Chr1SNPS = prettybaseread('1.prettybase.txt')
%

%   Copyright 2004 The MathWorks, Inc.

% More details on the Prettybase format are available from
% https://innateimmunity.net/Documentation/Prettybase
%
% Prettybase format SNPs for the human genome are available from
% http://pga.gs.washington.edu/data_download.html


try
    [snpOffset,subjectID,allele1,allele2] = textread(filename,'%s%s%s%s');
catch
    if(exist(filename) == 2)
        error('%s does not seem to be in PrettyBase format.',filename);
    else
        error('%s is not a valid file.',filename);
    end
end

if nargout < 2
    % output as a struct
    snpOffset = struct('snpOffset',snpOffset,'subjectID',subjectID,...
        'allele1',allele1,'allele2',allele2);
end

