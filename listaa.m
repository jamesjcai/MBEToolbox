function listaa

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


disp(' ')
disp('======================================================')
disp('Name           Abbr.          Linear structure formula')
disp('------------------------------------------------------')
disp(' ')
disp('Alanine        ala a                  CH3-CH(NH2)-COOH')
disp(' ')
disp('Arginine       arg r  HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH')
disp(' ')
disp('Asparagine     asn n           H2N-CO-CH2-CH(NH2)-COOH')
disp(' ')
disp('Aspartic acid  asp d             HOOC-CH2-CH(NH2)-COOH')
disp(' ')
disp('Cysteine       cys c               HS-CH2-CH(NH2)-COOH')
disp(' ')
disp('Glutamine      gln q        H2N-CO-(CH2)2-CH(NH2)-COOH')
disp(' ')
disp('Glutamic acid  glu e          HOOC-(CH2)2-CH(NH2)-COOH')
disp(' ')
disp('Glycine        gly g                      NH2-CH2-COOH')
disp(' ')
disp('Histidine      his h       NH-CH=N-CH=C-CH2-CH(NH2)-COOH')
disp('                           |__________|')
disp(' ')
disp('Isoleucine     ile i      CH3-CH2-CH(CH3)-CH(NH2)-COOH')
disp(' ')
disp('Leucine        leu l        (CH3)2-CH-CH2-CH(NH2)-COOH')
disp(' ')
disp('Lysine         lys k           H2N-(CH2)4-CH(NH2)-COOH')
disp(' ')
disp('Methionine     met m         CH3-S-(CH2)2-CH(NH2)-COOH')
disp(' ')
disp('Phenylalanine  phe f               Ph-CH2-CH(NH2)-COOH')
disp(' ')
disp('Proline        pro p                 NH-(CH2)3-CH-COOH')
disp('                                     |_________|')
disp(' ')
disp('Serine         ser s               HO-CH2-CH(NH2)-COOH')
disp(' ')
disp('Threonine      thr t           CH3-CH(OH)-CH(NH2)-COOH')
disp(' ')
disp('Tryptophan     trp w       Ph-NH-CH=C-CH2-CH(NH2)-COOH')
disp('                            |_______|')
disp(' ')
disp('Tyrosine       tyr y          HO-p-Ph-CH2-CH(NH2)-COOH')
disp(' ')
disp('Valine         val v            (CH3)2-CH-CH(NH2)-COOH')
disp('======================================================')
disp(' ')


% i_showimage('private/aminoacid.gif')

function i_showstereo(aa)
s=aa;
switch (lower(s))
    case ({'a','ala',1})
	disp('Alanine        ala a                  CH3-CH(NH2)-COOH')
	disp(' ')
	i_showimage('private/ala3d.gif')
end


function i_showimage(filename,formt)
	if (nargin<2), formt='gif'; end
	[X,MAP] = imread(filename,formt);
	colormap(MAP);
	image(X)
	daspect([1 1 1])
