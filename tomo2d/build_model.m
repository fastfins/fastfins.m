function model = build_model(n_cells, angle_sources, n_sources, n_detectors, detector_width)
%MODEL = BUILD_MODEL(n_cells, angle_sources, n_sources, n_detectors, detector_width)
%
% The model domain (object to be estimated) is defined in a box [-20,20]^2, 
% the radiation sources and detectors are located on a ring with radius 40
%
% n_cells:        number of discretization cells in 1D for the object
% angle_sources:  the span angle of radiation sources
% n_sources:      number of radiation sources
% n_detectors:    number of detectors per radiation source
% detector_width: width of each detector


model.B = 20; % box size
model.R = 40; % ring size

% fixed parameters
model.n_cells = 32;   % number of cells each direction
model.a_sourc = pi/2; % angle of radiation sources
model.n_sourc = 10;   % number of radiation sources
model.n_detec = 10;   % number of detectors
model.d_width = model.R*0.05; % detector width

switch nargin
    case {1}
        model.n_cells = n_cells;
    case {2}
        model.n_cells = n_cells;
        model.a_sourc= angle_sources;
    case {3}
        model.n_cells = n_cells;
        model.a_sourc = angle_sources;
        model.n_sourc = n_sources;
    case {4}
        model.n_cells = n_cells;
        model.a_sourc = angle_sources;
        model.n_sourc = n_sources;
        model.n_detec = n_detectors;
    case {5}
        model.n_cells = n_cells;
        model.a_sourc = angle_sources;
        model.n_sourc = n_sources;
        model.n_detec = n_detectors;
        model.d_width = detector_width;
end

% set the reference panel                  
an = atan(2*model.R/model.d_width) + pi/2;
th = 2*pi-2*an; % incremental angle
rot = [cos(th) -sin(th);sin(th) cos(th)];
panel = [cos(th/2); sin(th/2)]*model.R;
%
for i = 2:(model.n_detec/2)
    panel = [panel rot*panel(:,i-1)];
end
%
tmp = panel; 
tmp(2,:) = -tmp(2,:);
model.panel = [fliplr(panel) tmp];
% end of the reference panel 

if abs(model.a_sourc - 2*pi) < 1E-6
    as = linspace(-model.a_sourc/2, model.a_sourc/2, model.n_sourc+1);
    as = as(1:model.n_sourc);
elseif model.a_sourc < 2*pi
    as = linspace(-model.a_sourc/2, model.a_sourc/2, model.n_sourc);
else
    error('source angle should be less than or equal to 2*pi')
end
model.sources = [model.R*cos(as);model.R*sin(as)];
model.detectors = zeros(2, model.n_detec, model.n_sourc);
for i = 1:model.n_sourc
    % calculate detector angles, version 1
    % mid angle: as(i)+pi
    % start angle: mid angle - ang_detec
    % end angle: mid angle + ang_detec
    % ad = linspace(as(i)+pi-model.a_detec, as(i)+pi+model.a_detec, model.n_detec);
    % model.detectors(:,:,i) = [R*cos(ad);R*sin(ad)];
    % version 2
    model.detectors(:,:,i) = [cos(as(i)+pi) -sin(as(i)+pi); sin(as(i)+pi) cos(as(i)+pi)]*model.panel;
end

A = zeros(model.n_detec*model.n_sourc, model.n_cells^2);
%
grid1 = ((0:model.n_cells)/model.n_cells)*(2*model.B) - model.B;
for i = 1:model.n_sourc 
    P1 = model.sources(:,i);
    for j = 1:model.n_detec
        P2 = model.detectors(:,j,i);
        slope = (P1(2) - P2(2))/(P1(1) - P2(1)); 
        % intersection with x grid, gives x coord
        inter_x_grid = slope*(grid1 - P1(1)) + P1(2);
        % intersection with y grid, gives y coord
        inter_y_grid = (grid1 - P1(2))/slope + P1(1);
        % row lexicographical ordering 
        for row = 1:model.n_cells
            for col = 1:model.n_cells
                dist = cut_length(grid1(row), grid1(row+1), grid1(col), grid1(col+1), ...
                    inter_y_grid(row), inter_y_grid(row+1), inter_x_grid(col), inter_x_grid(col+1));
                %
                A(j+(i-1)*model.n_detec, col+(row-1)*model.n_cells) = dist;
            end
        end
    end
end

model.F = sparse(A);

figure('position', [100, 100, 1000, 500])
subplot(1,2,1)
hold on;
axis equal;
grid on
N = 1000;
as = linspace(0, 2*pi, N);
plot(model.R*cos(as), model.R*sin(as))
hold on
for i = 1:(model.n_cells+1)
    line([grid1(i),grid1(i)], [-model.B, model.B])
    line([-model.B, model.B], [grid1(i),grid1(i)])
end
i = 1;
plot(model.sources(1,i), model.sources(2,i), 'ro')
plot(model.detectors(1,:,i),model.detectors(2,:,i),'k.')
line([model.sources(1,i),model.detectors(1,1,i)], ...
    [model.sources(2,i), model.detectors(2,1,i)], 'color', 'r')
line([model.sources(1,i),model.detectors(1,end,i)], ...
    [model.sources(2,i), model.detectors(2,end,i)], 'color', 'r')
axis([-1,1,-1,1]*model.R)

subplot(1,2,2)
hold on;
axis equal;
grid on
N = 1000;
as = linspace(0, 2*pi, N);
plot(model.R*cos(as), model.R*sin(as))
hold on
for i = 1:(model.n_cells+1)
    line([grid1(i),grid1(i)], [-model.B, model.B])
    line([-model.B, model.B], [grid1(i),grid1(i)])
end
plot(model.sources(1,:), model.sources(2,:), 'ro')
for i = 1:model.n_sourc
    plot(model.detectors(1,:,i),model.detectors(2,:,i),'k.')
    line([model.sources(1,i),model.detectors(1,1,i)], ...
        [model.sources(2,i), model.detectors(2,1,i)], 'color', 'r')
    line([model.sources(1,i),model.detectors(1,end,i)], ...
        [model.sources(2,i), model.detectors(2,end,i)], 'color', 'r')
end
axis([-1,1,-1,1]*model.R)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist = cut_length(in, i1n, jn, j1n, xin, xi1n, yjn, yj1n)

% from the old code

d = 0;
if (in < yjn) && (yjn < i1n) 
    if (jn < xin) && (xin < j1n) 
        d = sqrt((jn - xin)^2 + (yjn - in)^2);
    end
    if (in < yj1n) && (yj1n < i1n)  
        d = sqrt((jn - j1n)^2 + (yjn - yj1n)^2);
    end
    if (jn < xi1n) && (xi1n < j1n) 
        d = sqrt((jn - xi1n)^2 + (yjn - i1n)^2);
    end
end
if (jn < xin) && (xin < j1n)    
    if (jn < xi1n) && (xi1n < j1n) 
       d = sqrt((xin - xi1n)^2 + (in-i1n)^2);
    end
    if (in < yj1n) && (yj1n < i1n)  
       d = sqrt((xin - j1n)^2 + (in - yj1n)^2);    
    end
end
if (jn < xi1n) && (xi1n < j1n) 
    if (in < yj1n) && (yj1n < i1n) 
        d = sqrt((xi1n - j1n)^2 + (i1n - yj1n)^2);
    end
end                        


if (yjn == i1n) && (jn < xin) && (xin < j1n) 
    d = sqrt((jn - xin)^2 + (in - i1n)^2);
end
if (yjn == i1n) && (yj1n == in)  
    d = sqrt((jn - j1n)^2 + (in - i1n)^2);
end
if (yjn == i1n) && (in < yj1n) && (yj1n < i1n) 
    d = sqrt((jn - j1n)^2 + (i1n - yj1n)^2);
end
if (yjn == in) && (jn < xi1n) && (xi1n < j1n) 
    d = sqrt((jn - xi1n)^2 + (in - i1n)^2);
end
if (yjn == in) && (in < yj1n) && (yj1n < i1n) 
    d = sqrt((jn - j1n)^2 + (in - yj1n)^2);
end
if (yjn == in) && (yj1n == i1n)   
    d = sqrt((in - i1n)^2 + (jn - j1n)^2);
end
if (yj1n == in) && (in < yjn) && (yjn < i1n)  
    d = sqrt((jn - j1n)^2 + (in - yjn)^2);
end
if (yj1n == in) && (jn < xi1n) && (xi1n < j1n) 
    d = sqrt((in - i1n)^2 + (j1n - xi1n)^2);
end
if (yj1n == i1n) && (in < yjn) && (yjn < i1n)  
    d = sqrt((jn - j1n)^2 + (i1n - yjn)^2);
end
if (yj1n == i1n) && (jn < xin) && (xin < j1n)
    d = sqrt((i1n - in)^2 + (j1n - xin)^2);
end

dist = d;

end