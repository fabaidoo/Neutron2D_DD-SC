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
        
        function [psi_out, Qnew] = stepchar(obj, Ox, psi_in, phi0)
            %takes incoming angular flux and previous 0th and 1st angular 
            %moments of average and calculates outgoing angular flux via
            %diamond difference method
            Delta = obj.sidelength;
            tau = obj.sig_t * Delta ./ abs(Ox); 
            
            %total source term
            Qnew = obj.Q + 0.5 * obj.sig_s0 * phi0;
            
            %compute the exiting angular flux
            psi_out = psi_in .* exp(-tau) + (Qnew ./ obj.sig_t) .* ...
                (1 - exp(-tau)); 
        end
        
        
         function [psi_out, Qnew] = diamdiff(obj, Ox, psi_in, phi0)
            %takes incoming angular flux and previous 0th and 1st angular 
            %moments of average and calculates outgoing angular flux via
            %diamond difference method
            Delta = obj.sidelength;
            tau = obj.sig_t * Delta ./ abs(Ox); 
            
            %total source term
            Qnew = obj.Q + 0.5 * obj.sig_s0 * phi0;
            
             %compute the exiting angular flux
            psi_out = psi_in .* (2 - tau) ./ (2 + tau) + ... 
                (Qnew ./ obj.sig_t) .* (1 - (2 - tau) ./ (2 + tau)); 
         end
        
         function phi0_out = phi_maker(obj, Ox, w, Q, psi0,...
                psi1)
            %Provides terms for new angular moments of phi for material. 
            %Takes angle Oz and its weight, incoming and outgoing psis
            %and modified source
            
            Delta = obj.sidelength; %material width
            tau = obj.sig_t * Delta ./ abs(Ox); %has same shape as Oz 
            
            phi0_out =  (Q / obj.sig_t + (psi0 - psi1) / tau ) * w;
            
           
        end
            
         
         
         
    end
end

