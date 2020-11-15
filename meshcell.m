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
                obj.sig_s = 2;
            elseif strcmpi(obj.material, 'air') == 1
                obj.sig_t = 0.01;
                obj.sig_s = 0.006;
            else
                error('material type not available')
            end
            
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
            
            f = @(x) (exp(- x) - 1) / x;
            
            psi0_ave = S * (1 /rho + rho) / 2 + (psix_in - S) * (1 + f(a)) /a + (psiy_in - S) * (1 + f(b)) /b; 
            
        end
  
    end
    
    methods(Static)
         function psi_out = diamdiff(psi_in, psi_ave)
            %takes incoming angular flux and average angular flux and 
            %mcalculates outgoing angular flux for diamond difference 
            %method
            psi_out = 2 * psi_ave - psi_in; 
         end
      
    end  
end

