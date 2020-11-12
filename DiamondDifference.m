function DiamondDifference(meshnum, n_xy, n_z)

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

[Ox, Oy, w] = angle_x(n_xy, n_z); %angular discretization

%initial guesses for moment of scalar flux
phi0 = zeros(meshnum);

%angular fluxes go here (VACUUM BCs)
psix_l = zeros(meshnum, meshnum + 1) ;%forward sweeps
psix_r = zeros(meshnum, meshnum + 1) ;%backward sweeps

psiy_b = zeros(meshnum + 1, meshnum); %forward sweeps
psiy_t = zeros(meshnum + 1, meshnum); %backward sweeps


max_iter = 3e02; %break while loop after max_iter iterations 
iter = 0; %iteration count
err = .0050  ; %error in calculated phis
while iter <= max_iter && err > tol
    phi0old = phi0 ;
     
    for i = 1:length(Ox)
        if Ox(i) > 0
            for j = 1: length(Oy)
                if Oy(j) > 0
                    for k = 1: meshnum 
                        for l = 1: meshnum 
                            obj = mesh{k,l};
                            [psix_l(k, l+1), ~] = obj.diamdiff(Ox(i), psix_l(k,l), phi0(k,l));
                            [psiy_b(k+1, l), Q] = obj.diamdiff(Oy(j), psiy_b(k, l), phi0(k,l));
                            
                            phi0(k,l) = obj.phi_maker(Ox(i), w, Q, psix_l(k,l), psix_l(k,l+1)); %HAVE TO UPDATE PHI_MAKER METHOD!!!!!
                        end
                    end
                end
                
                if Oy(j) < 0
                    for k = 1: meshnum
                        
                        for l = 1: meshnum 
                            m = meshnum + 1 - k; %!!!!!!!!!!!
                            obj = mesh{m,l};
                            [psix_l(m, l+1), ~] = obj.diamdiff(Ox(i), psix_l(m,l), phi0(m,l));
                            [psiy_t(m, l), Q] = obj.diamdiff(Oy(j), psiy_t(m+1, l), phi0(m,l));
                            
                            phi0(m,l) = obj.phi_maker(Ox(i), w, Q, psix_l(m,l), psix_l(m,l+1)); %HAVE TO UPDATE PHI_MAKER METHOD!!!!!
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
                            [psix_r(k, n), ~] = obj.diamdiff(Ox(i), psix_r(k,n+1), phi0(k,n));
                            [psiy_b(k+1, n), Q] = obj.diamdiff(Oy(j), psiy_b(k, n), phi0(k,n));
                            
                            phi0(k,n) = obj.phi_maker(Ox(i), w, Q, psix_l(k,n+1), psix_l(k,n)); %HAVE TO UPDATE PHI_MAKER METHOD!!!!!
                        end
                    end
                end
                
                if Oy(j) < 0
                    for k = 1: meshnum
                         for l = 1: meshnum
                             m = meshnum + 1 - k;
                             n = meshnum + 1 - l;
                             obj = mesh{m,n};
                             [psix_r(m, n), ~] = obj.diamdiff(Ox(i), psix_r(m,n+1), phi0(m,n));
                             [psiy_t(m, n), Q] = obj.diamdiff(Oy(j), psiy_t(m+1, n), phi0(m,n));
                            
                             phi0(m,n) = obj.phi_maker(Ox(i), w, Q, psix_l(m,n+1), psix_l(m,n)); %HAVE TO UPDATE PHI_MAKER METHOD!!!!!
                         end
                    end
                end
            end
        
        end
    end
   err = norm(phi0old - phi0norm, 'fro');
   iter = iter + 1;
    
end


%{
phi_diag = diag(phi0);
phi_bottom = phi0(1, :);

figure
plot(cent, phi_bottom)

figure
plot(cent, phi_diag)
%}


end