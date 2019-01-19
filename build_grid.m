function[Grid]=build_grid(Grid)
%author: Peter Carlson
%date: 2/3/2015
%description:
%This function accepts a minimal definition of the computational domain
%and grid and computes all pertinent information about the grid.
%Input:
%Grid.xmin = left boundary of the domain
%Grid.xmax = right boundary of the domain
%Grid.Nx = number of grid cells
%Output:
%Grid.Lx = Length of the domain
%Grid.dx = cell width
%Grid.xc = vector of cell center locations
%Grid.xf = vector of cell face locations
%Grid.Nfx = number of fluxes in x-directions
%Grid.dof = vector from 1 to N containing degrees of freedom
%Grid.dof_xmin = degrees of freedom corresponding to the left boundary
%Grid.dof_xmax = degrees of freedom corresponding to the right boundary
%Grid.geom = coordinate system to use: cartesian = 1, polar = 2, 
%spherical = 3
%
%Example Call:
%>>Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 5;
%>>Grid = build_grid(Grid);


Grid.Lx = Grid.xmax - Grid.xmin;
Grid.dx = Grid.Lx/Grid.Nx;
Grid.xc = (Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax-Grid.dx/2);
Grid.xf = (Grid.xmin:Grid.dx:Grid.xmax);
Grid.Nfx = Grid.Nx+1;
Grid.dof = Grid.Nx;
Grid.dof_xmin = Grid.xmin;
Grid.dof_xmax = Grid.Nx;

isfieldresult = isfield(Grid, 'geom');
if isfieldresult == 0
    Grid.geom = 1;
end

