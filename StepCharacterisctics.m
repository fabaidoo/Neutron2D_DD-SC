function StepCharacterisctics(meshnum)

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




end