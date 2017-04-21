%% Converts a 3-dimensional RGB image from RGB to LAB colorspace.
function [ A ] = colorspace3d( I, rgb2lab )

if nargin < 2
    rgb2lab = 1;
end

dim = size(I);
A = zeros(dim);

    if rgb2lab
        for s = 1:dim(3)
            J = squeeze(I(:,:,s,:));
            % Convert input sRGB image to CIELab color space.
            A(:,:,s,:) = colorspace('Lab<-RGB',J);
        end
    else
        for s = 1:dim(3)
            J = squeeze(I(:,:,s,:));
            % Convert filtered image back to sRGB color space.
            A(:,:,s,:) = colorspace('RGB<-Lab',J);
        end
    end
end

