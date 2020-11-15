function info = DiamondDifference(meshnum, n_xy, n_z)
innerbox = .1; % size of inner box with source
outerbox = 1; %size of outer box
material = 'absorber'; %material in outer box 


%CREATE DOMAIN AND MESH
h = outerbox/meshnum; %meshsize
edges = 0: h : outerbox; %ALSO: linspace(0,outerbox, meshnum + 1)

cent = edges(1:meshnum) + diff(edges)/2; %coordinates of mesh elements

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

[Ox, Oy, w] = angle(n_xy, n_z); %angular discretization

%initial guesses for moment of scalar flux
phi0 = zeros(meshnum);
%phi0old = rand(meshnum);
%angular fluxes go here (VACUUM BCs)
psix_l = zeros(length(Ox), meshnum, meshnum + 1) ;%forward sweeps
psix_r = zeros(length(Ox), meshnum, meshnum + 1) ;%backward sweeps

psiy_b = zeros(length(Ox), meshnum + 1, meshnum); %forward sweeps
psiy_t = zeros(length(Ox), meshnum + 1, meshnum); %backward sweeps


max_iter = 3e02; %break while loop after max_iter iterations 
iter = 0; %iteration count
err = 50  ; %error in calculated phis
tol = 1e-5;

phi0_ang = zeros(length(Ox), meshnum, meshnum); %cell average psi goes here

while iter <= max_iter && err > tol
    
    for i = 1:length(Ox)
        if Ox(i) > 0
            if Oy(i) > 0
                for k = 1: meshnum
                    for l = 1: meshnum
                        obj = mesh{k,l};
                        [psi0_ave, psix_l(i, k, l+1), psiy_b(i, k+1, l) ] =...
                            obj.diamonddifference(Ox(i), Oy(i),  psix_l(i, k,l), psiy_b(i, k, l), phi0(k,l));
                        
                        phi0_ang(i, k,l) = psi0_ave * w(i);
                    end 
                end
            end
                
            if Oy(i) < 0
                for k = 1: meshnum
                        
                    for l = 1: meshnum
                        m = meshnum + 1 - k;
                        obj = mesh{m,l};
                        [psi0_ave, psix_l(i, m, l+1), psiy_t(i, m, l)]  = ...
                            obj.diamonddifference(Ox(i), Oy(i),  psix_l(i, m,l), psiy_t(i, m+1, l), phi0(m,l));
                        
                        phi0_ang(i, m,l) = psi0_ave * w(i);
                    end
                end
            end
        end
        
        if Ox(i) < 0
            if Oy(i) > 0
                for k = 1: meshnum
                    for l = 1: meshnum
                        n = meshnum + 1 - l;
                        obj = mesh{k,n};
                        [psi0_ave, psix_r(i, k, n), psiy_b(i, k+1, n)] = ...
                            obj.diamonddifference(Ox(i), Oy(i), psix_r(i, k,n+1), psiy_b(i, k, n), phi0(k,n));
                        
                        phi0_ang(i, k, n) = psi0_ave * w(i);
                    end
                end
            end
                
            if Oy(i) < 0
                for k = 1: meshnum
                    for l = 1: meshnum
                        m = meshnum + 1 - k;
                        n = meshnum + 1 - l;
                        obj = mesh{m,n};
                        [psi0_ave, psix_r(i, m, n), psiy_t(i, m, n)] = ...
                            obj.diamonddifference(Ox(i), Oy(i), psix_r(i, m,n+1), psiy_t(i, m+1, n), phi0(m,n));
                        
                        phi0_ang(i, m,n) = psi0_ave * w(i);
                    end
                end
                
            end
        
        end
    end
    
    phi0old = phi0;
    phi0 = reshape(sum(phi0_ang), [meshnum, meshnum]);
    err = norm(phi0old - phi0, 'fro');
    iter = iter + 1;
   
  % 
   if iter == max_iter
       error('Maximum number of iterations reached before convergence. Error = %.6f', err)
   end
   
    if err > 5e07 
        disp(err)
        error('Solution is blowing up. iter = %i', iter)
      
   end
  %}  
end

phi_diag = diag(phi0);
phi_bottom = phi0(:, end);

str = sprintf('%s | Space: %i cells/edge. | Angular: m_{xy} = %i, m_z = %i ', upper(material), meshnum,  n_xy, n_z);

figure
plot(cent, phi_bottom, '.-r','MarkerSize', 11, 'LineWidth', 1.5)
xlabel('(x, 0)')
ylabel('\phi(x, 0)')
title(str)

figure
plot(cent, phi_diag, '.-', 'MarkerSize', 11, 'LineWidth', 1.5)
ylabel('\phi(x, x)')
xlabel('(x, x)')
title(str)

figure
surf(cent, cent, phi0)
xlabel('x')
ylabel('y')
zlabel('\phi(x,y)')
title(str)

%}

info = sprintf("iter = %i \nerror = %.3e ", iter, err);




end