% ChenVarshney.m
% -------------------------------------------------------------------
% Reference: <A human perception inspired quality metric for image 
%               fusion based on regional information>
%            <Objective assessment of multiresolution image fusion algorithms 
%               for context enhancement in Night vision: A comparative study>
%             http://www.mathworks.com/matlabcentral/fileexchange/3689-wpsnr/content/wpsnr.m
%             http://vision.arc.nasa.gov/personnel/al/code/matlab/filtmod2.htm
% Authors: Sun Li
% Date:    03/12/2014
% Last modified: 03/12/2014
% -------------------------------------------------------------------


function F = CSFFilterF(img, R)

    ff = fft2(img);
    ffc = fftshift(ff);
    [m, n] = size(img);
    [xx, yy] = meshgrid(1:n, 1:m);
    xx = xx - (n+1)/2;
    yy = yy - (m+1)/2;
    
    xx = xx/n*R;
    yy = yy/m*R;
    fw = MannosSkarision(sqrt(xx.^2+yy.^2));
    %fw = ones(m,n);
    affc = abs(ffc(:));
    energ = sqrt(sum(affc.*affc));
    %ffc = ffc/(asum+eps);
    %b = sum(sum((abs(ffc)).^2));
    F = abs(ffc.*fw);
    %F = sum(F(:))/energ;
    F = sqrt(sum(F(:).*F(:)))/energ;
    %F = sqrt(sum(F(:).*F(:))/(m*n));
    %F = sum(F(:).*F(:))/energ;
end


%%
function fc = MannosSkarision(r)     
    fc = 2.6*(0.0192+0.114*r).*exp(-(0.114*r).^1.1); % The same as paper
end


function fc = AhumadaFilter(r)
    ac = 1;
    as = 0.685;
    fc = 97.3227;
    fs = 12.1653;
    fc = ac*exp(-(r/fc).^2)-as*exp(-(r/fs).^2);
end