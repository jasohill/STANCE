function CenterOfMass = center_of_mass(data)
% Find the center of mass coordinates of a 3d data array.

   size_data = size(data);
   mass = sum(data(:));
   CenterOfMass = zeros(3,1);

   for di = 1:3
       size_di = ones(1,3);
       size_di(di) = size_data(di);
       replicate = size_data;
       replicate(di) = 1;
       indices = repmat(reshape(1:size_data(di),size_di),replicate);
       CenterOfMass(di) = sum(indices(:).*data(:))./mass;
   end

end

