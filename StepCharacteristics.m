function info = StepCharacteristics(meshnum, n_xy, n_z)

innerbox = .1; % size of inner box with source
outerbox = 1; %size of outer box
material = 'absorber'; %material in outer box 


%create domain and mesh
h = outerbox/meshnum; %meshsize
edges = 0: h : outerbox; %linspace(0,outerbox, meshnum + 1)
%coordinates of center of mesh elements
cent = edges(1:meshnum) + diff(edges)/2; 

mesh = cell(meshnum); %mesh elements go here
for i = 1:meshnum
    for j = 1: meshnum
        center = [cent(i), cent(j)];
        x_dist = abs(cent(i) - outerbox/2);
        y_dist = abs(cent(j) - outerbox/2);
        if max(x_dist, y_dist) < innerbox/2
            mat = 'source';
        else
            mat = material;
        end
        
        mesh{i,j} = meshcell(mat, center, h);     
    end 
end

[Ox, Oy, w] = angle_x(n_xy, n_z); %angular discretization

%initial guesses for moment of scalar flux
phi0 = zeros(meshnum);
%phi0old = rand(meshnum);
%angular fluxes go here (VACUUM BCs)
psix_l = zeros(length(Ox), length(Oy), meshnum, meshnum + 1) ;%forward sweeps
psix_r = zeros(length(Ox), length(Oy), meshnum, meshnum + 1) ;%backward sweeps

psiy_b = zeros(length(Ox), length(Oy), meshnum + 1, meshnum); %forward sweeps
psiy_t = zeros(length(Ox), length(Oy), meshnum + 1, meshnum); %backward sweeps


max_iter = 3e02; %break while loop after max_iter iterations 
iter = 0; %iteration count
err = 50  ; %error in calculated phis
tol = 1e-5;

phi0_ang = zeros(length(Ox), length(Oy), meshnum, meshnum);

while iter <= max_iter && err > tol
    for i = 1:length(Ox)
        if Ox(i) > 0
            for j = 1: length(Oy)
                if Oy(j) > 0
                    for k = 1: meshnum 
                        for l = 1: meshnum 
                            obj = mesh{k,l};
                            [psi0_ave, psix_l(i, j, k, l+1), psiy_b(i, j, k+1, l)] = obj.stepcharacteristics(Ox(i), Oy(j), ...
                                psix_l(i, j, k,l), psiy_b(i, j, k, l), phi0(k,l));
                            phi0_ang(i,j,k,l) = psi0_ave; %phi0(k,l) + psi0_ave;
                            
                        end
                    end
                end
                
                if Oy(j) < 0
                    for k = 1: meshnum
                        
                        for l = 1: meshnum 
                            m = meshnum + 1 - k; %!!!!!!!!!!!
                            obj = mesh{m,l};
                            [psi0_ave,psix_l(i, j, m, l+1), psiy_t(i, j, m, l)] = obj.stepcharacteristics...
                                (Ox(i), Oy(j),  psix_l(i, j, m,l), psiy_t(i, j, m+1, l), phi0(m,l));
                            phi0_ang(i,j, m,l) = psi0_ave;   %phi0(m,l) + psi0_ave;
                        end
                    end
                end
            end
        end
        
        if Ox(i) < 0
            for j = 1: length(Oy)
                if Oy(j) > 0
                    for k = 1: meshnum
                        for l = 1: meshnum 
                            n = meshnum + 1 - l;
                            obj = mesh{k,n};
                            [psi0_ave,psix_r(i, j, k, n),psiy_b(i, j, k+1, n)] = obj.stepcharacteristics...
                                (Ox(i), Oy(j), psix_r(i, j, k,n+1),psiy_b(i, j, k, n), phi0(k,n));
                            phi0_ang(i,j, k, n) = psi0_ave;%phi0(k,n) + psi0_ave;
                        end
                    end
                end
                
                if Oy(j) < 0
                    for k = 1: meshnum
                         for l = 1: meshnum
                             m = meshnum + 1 - k;
                             n = meshnum + 1 - l;
                             obj = mesh{m,n};
                             [psi0_ave, psix_r(i, j, m, n), psiy_t(i, j, m, n)] = obj.stepcharacteristics...
                                 (Ox(i), Oy(j), psix_r(i, j, m,n+1), psiy_t(i, j, m+1, n), phi0(m,n));
                             phi0_ang(i,j, m,n) = psi0_ave;% phi0(m,n) + psi0_ave;
                         end
                    end
                end
            end
        
        end
    end
   phi0old = phi0;
   phi0 = reshape(w * sum(sum(phi0_ang)), [meshnum, meshnum]); 
   err = norm(phi0old - phi0, 'fro');
   iter = iter + 1;
   
  % 
   if iter == max_iter
       error('Maximum number of iterations reached before convergence. Error = %.6f', err)
   end
   
   if err > 5e030
       disp(err)
       error('Solution is blowing up. iter = %i', iter)
   end
  %}  
end


%
%phi_diag = diag(phi0);
%phi_bottom = phi0(:, end);

%figure
%plot(cent, phi_bottom)

%figure
%plot(cent, phi_diag, '*')

figure
surf(cent, cent, phi0)

%}


info = sprintf("iter = %i \nerror = %.3g ", iter, err);




end