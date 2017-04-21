function M2T = mni2tal_matrix()
% mni2tal_matrix - Talairach to MNI coordinates (best guess)
%
% MNI2TALAIRACH = mni2tal_matrix
%
% MNI2TALAIRACH is a struct containing rotation matrices
% used by mni2tal and tal2mni
%
% See also, MNI2TAL, TAL2MNI &
% http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
%
% $Revision: 1.2 $ $Date: 2004/11/17 21:04:32 $
% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 2/2/01, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - removed dependence on spm_matrix by
%                     creating this function, thereby
%                     abstracting the important matrix
%                     transforms (easier to change if needed).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See notes below for explanations...

% rotn  = spm_matrix([0 0 0 0.05]); % similar to Rx(eye(3),-0.05), DLW
M2T.rotn  = [      1         0         0         0;
                   0    0.9988    0.0500         0;
                   0   -0.0500    0.9988         0;
                   0         0         0    1.0000 ];

% upz   = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.92]);
M2T.upZ   = [ 0.9900         0         0         0;
                   0    0.9700         0         0;
                   0         0    0.9200         0;
                   0         0         0    1.0000 ];

% downz = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.84]);
M2T.downZ = [ 0.9900         0         0         0;
                   0    0.9700         0         0;
                   0         0    0.8400         0;
                   0         0         0    1.0000 ];

% from original mni2tal...
%upT   = spm_matrix([0 0 0 0.05 0 0 0.99 0.97 0.92]);
%downT = spm_matrix([0 0 0 0.05 0 0 0.99 0.97 0.84]); 

return