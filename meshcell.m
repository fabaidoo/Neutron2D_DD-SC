classdef meshcell
    %CONSTITUENT MESH ELEMENT 
    % Must be made of one type of material. Class holds material 
    %properties as well as location in space.
    
    properties
        sig_t = 0; %TOTAL CROSS-SECTION
        sig_s = 0; %SCATTERING CROSS-SECTION (0TH MOMENT)
        Q = 0; %SOURCE TERM (0TH MOMENT)
        
        center = [0 0]; %LOCATION OF CENTER OF SQUARE ELEMENT
        sidelength = 1; %LENGTH OF SQUARE ELEMENT
        
        material = 'source'; %TYPE OF MATERIAL. Helpful for debugging
    end
    
    methods
        function obj = meshcell(material,center, sidelength)
            %Specify type of material and location in space.
            if nargin == 1
                obj.material = material;
            elseif nargin == 3
                obj.material = material;
                obj.center = center;
                obj.sidelength = sidelength;
            end
            
            if strcmpi(obj.material, 'source') == 1
                obj.sig_t = 0.1;
                obj.Q = 1;
            elseif strcmpi(obj.material, 'scatterer') == 1
                obj.sig_t = 2;
                obj.sig_s = 1.99;
            elseif strcmpi(obj.material, 'reflector') == 1
                obj.sig_t = 2;
                obj.sig_s = 1.8;
            elseif strcmpi(obj.material, 'absorber') == 1
                obj.sig_t = 10;
                obj.sig_s = 0;%2; %!!!!!!!!!
            elseif strcmpi(obj.material, 'air') == 1
                obj.sig_t = 0.01;
                obj.sig_s = 0.006;
            else
                error('material type not available')
            end
            
        end
        
        function  [psi_ave, psix_out, psiy_out]= stepchar(obj, Ox, Oy, psix_in, psiy_in, phi0old)
            %takes incoming angular flux and previous 0th and 1st angular 
            %moments of average and calculates outgoing angular flux via
            %diamond difference method
            h = obj.sidelength;
            h_real = [h/abs(Ox), h/abs(Oy)];
            [d, p] = min(h_real); %distance and index telling us if min is from x or y
            q = 3 - p; %index of the max term above 
            
            %total source term
            Qnew = obj.Q + 0.5 * obj.sig_s * phi0old;
            
            psi_in = [psix_in, psiy_in]; %incoming psi
            
            psi_0 = Qnew/ obj.sig_t + obj.p_exp(-obj.sig_t * d, 4) * (psi_in - Qnew / obj.sig_t);
            psi_1 = Qnew / obj.sig_t + (psi_in - psi_0) / (obj.sig_t * d);
            
            %normalized distance
            del = zeros(1, 2);
            del(p) = d / h_real(p);
            del(q) = d / h_real(q);
            
            psi_out = zeros(1, 2);
            psi_out(p) = del(q) * psi_1(q) + psi_0(p);
            psi_out(q) = del(p) * psi_1(p);
            
            psi_ave = Qnew / obj.sig_t + ((psi_in(p) - psi_out(p)) + del(q) * (psi_in(q) - psi_out(q)) ) / (d * obj.sig_t);
            
            psix_out = psi_out(1);
            psiy_out = psi_out(2);   
        end
      
        function psi0_ave = DDphi_maker(obj, Ox, Oy, psix_in,...
                 psiy_in, phi0_old)
            %Calculates average angular flux in cell for given x and y   
            %directions
            h = obj.sidelength; %material width
            Qnew = obj.Q + 0.5 * obj.sig_s * phi0_old; %total source term
            taux = 2 * abs(Ox) / h ;
            tauy = 2 * abs(Oy) / h ;
            
            psi0_ave = (Qnew + taux * psix_in + tauy * psiy_in) / ...
                (obj.sig_t + taux + tauy);
           
        end
        
        function [psi0_ave, psix_out, psiy_out] = diamonddifference(obj, Ox, Oy, psix_in, psiy_in, phi0_old)
             %Produces average psi in a cell and psi flow out of x and y
             %directions using diamond difference method
             
             h = obj.sidelength; %material width
             Qnew = obj.Q + 0.5 *  obj.sig_s * phi0_old; %total source term 
             taux = 2 * abs(Ox) / h ;
             tauy = 2 * abs(Oy) / h ;
            
             psi0_ave = (Qnew + taux * psix_in + tauy * psiy_in) / ...
                 (obj.sig_t + taux + tauy);
             
             psix_out = obj.diamdiff(psix_in, psi0_ave);
             psiy_out = obj.diamdiff(psiy_in, psi0_ave);
             
        end
        
        function [psi0_ave, psix_out, psiy_out] = stepcharacteristics(obj, Ox, Oy, psix_in, psiy_in, phi0_old)
            %Calculates outgoing fluxes and average flux in cell using step
            %characteristics method
            
            h = obj.sidelength;
            a = obj.sig_t * h / abs(Ox);
            b = obj.sig_t * h / abs(Oy);
            rho = a / b ;
            
            Qnew = obj.Q + 0.5 *  obj.sig_s * phi0_old;
            S = Qnew / obj.sig_t;
            
            if rho < 1
                psix_out = S + (psix_in - S) * (1 - rho) * exp(-a) +...
                    rho * (psiy_in - S) * (1 - exp(-a)) / a;
                
                psiy_out = S + (psix_in - S) * (1 - exp(-a)) / a ; 
                
            elseif rho >= 1
                psix_out = S + (psiy_in - S) * (1 - exp(-b) ) / b ;
                
                psiy_out = S + (psix_in - S) * (1 - exp(-b))/ (rho * b)+ ...
                    (psiy_in - S) * (1 - 1 /rho) * exp(-b) ; 
            end
            
            psi0_ave = S + (psix_out - psix_in) /a + (psiy_out - psiy_in) /b; 
        end
  
         
         function y = p_exp(obj, x, n)
             %taylor polynomial approximation to exp(x)
             y = 1;
             for i = 1:n
                y = y + obj.n_th_term(x, i); 
             end
         end
    end
    
    methods(Static)
         function psi_out = diamdiff(psi_in, psi_ave)
            %takes incoming angular flux and average angular flux and 
            %mcalculates outgoing angular flux for diamond difference 
            %method
            psi_out = 2 * psi_ave - psi_in; 
         end
         
         function y = n_th_term(x, n)
             %n-th term in taylor expansion of exp(x)
             y = x.^n ./ factorial(n);
             
         end
   
         
         
    end
    
    
    
    
end

