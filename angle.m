function [Ox, Oy, w] = angle(n_xy, n_z)
% xy - direction discretization along with weights 

if mod(n_xy, 4) ~= 0
    error('Quarter symmetry is required for the problem')
end
if mod(n_z, 2) ~= 0
   error('Only even discretizations in z are allowed') 
end

w_xy = 2 * pi / n_xy;
w_xy = w_xy *ones(1, n_xy);


phi = pi * linspace(1/n_xy , (2*n_xy - 1)/n_xy, n_xy);%azimuthal angle 
% SAME AS: 1/n: 2/n: (2*n - 1)/n

[cos_theta, w_z] = angle_z(n_z); 
sin_theta = sqrt(1 - cos_theta.^2 ); %polar component

Ox =  kron(cos(phi), sin_theta); %x-directions 
Oy = kron(sin(phi), sin_theta);

w = kron(w_xy, w_z);


function [mu, w_z] = angle_z(m)
%Angular discretization in the z-direction
if mod(m,2) == 1
    error('Even angular discretization needed')   
end
S2 = .5773502691896257;
w2 = 1;
       
S4 = [.3399810435848563 .8611363115940526];
w4 = [0.6521451548625461 0.3478548451374538];
      
S6 =[ .2386191860831969 .6612093864662645 .9324695142031521];
w6 = [0.4679139345726910 0.3607615730481386	 0.1713244923791704];
     
S8 = [.1834346 .5255324 .7966665 .9602899] ;
w8 = [.3626838 .3137066 .2223810 .1012285];
       
Sn = {S2, S4, S6, S8};
wn = {w2, w4, w6, w8};


mu = Sn{m/2};
w_z = wn{m/2};
end




end