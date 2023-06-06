%% RAND_SAMPLEMASK4K_ORI : generate  random sampling mask
%                   programmed by Wan Kim, May. 1, 2014
%
% randMask :    random mask for random sampling
% csrfac_real : real cs reduction factor measured from randMask

function rand2DMask = rand_sampleMask4k_ori(Ny,Nx,csrfac,alpha,radiu)

% Set parameters
mu = 0;      % the inverse of slope
delta = 1;   % increment of mu
iter = 1;    % iteration flag
VIP_dist = 1;

rand2Dnoise = rand(Ny,Nx);
rand2DMask = zeros(size(rand2Dnoise));
while(iter)
   mu = mu+delta;
   
   d0 = sqrt(((abs(Ny/2)).^2)+(abs(Nx/2)).^2); 
   d = sqrt(((abs([1:Ny]-Ny/2)).^2)'*ones(1,Nx)+ones(Ny,1)*(abs([1:Nx]-Nx/2)).^2); 

   threshold = exp(-1/mu*(d/d0).^alpha); % for Kim

   rand2DMask(rand2Dnoise<=threshold)=1;
   rand2DMask(rand2Dnoise>threshold)=0;
   rand2DMask(d<=VIP_dist)=1;
   
     N = size(rand2Dnoise);
     R = radiu;
    X_I = [floor(N(1)/2-R+1):floor(N(1)/2+R)];
    Y_I = [floor(N(2)/2-R+1):floor(N(2)/2+R)];
    vec_mask(X_I,Y_I) = ones(length(X_I),length(Y_I));
    [Ic,Jc,Sc] = find(vec_mask);
    iter = 1;
 while(iter <= length(Sc))
     D = sqrt((Ic(iter)-N(1)/2).^2+(Jc(iter)-N(2)/2).^2);
     if( D < R)
         rand2DMask(Ic(iter),Jc(iter)) = 1;
     end
     iter = iter + 1;
 end
   
   csrfac_real=1/mean(rand2DMask(:));
   err=csrfac_real-csrfac;
   if(err<0)
      mu=mu-delta;
      delta=delta/2;
      if(delta<1e-4)
       iter=0;
      end
   end
end

% figure;
% imshow(rand2DMask,[]);
% hold on; title('random pattern in k domain');
% h=sum(rand2DMask,2);
% figure; plot([0:Ny-1],h);
% xlabel('k-axis'); ylabel('histogram of 1s');
% axis([1, Ny,0, Nx]);
% 
