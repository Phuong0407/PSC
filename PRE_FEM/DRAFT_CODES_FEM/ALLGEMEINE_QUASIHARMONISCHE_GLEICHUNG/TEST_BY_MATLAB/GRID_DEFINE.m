% Define the polygon vertices in counterclockwise order
x = [0, 3, 6, 9, 9, 0];  % z-coordinates (x-axis in axisymmetric terms)
y = [0, 0, 2, 2, 4, 4];  % r-coordinates (y-axis in axisymmetric terms)

% Create the geometry description matrix
gd = [2; 6; x(:); y(:)];  % 2D polygon geometry
sf = 'P1';                % Set formula for single polygon
ns = char('P1');          % Name of the geometry
ns = ns';

% Decompose the geometry
dl = decsg(gd, sf, ns);

% Create an axisymmetric thermal model
model = createpde('thermal', 'steadystate-axisymmetric');

% Assign the decomposed geometry to the PDE model
geometryFromEdges(model, dl);

% Apply boundary conditions using thermalBC
% Edge 1: Lower boundary (z = 0, r = 0 → 3): u = 0
thermalBC(model, 'Edge', 1, 'Temperature', 0);

% Edge 2: Lower boundary (z = 3 → 6): u = 0
thermalBC(model, 'Edge', 2, 'Temperature', 0);

% Edge 3: Lower boundary (z = 6 → 9): u = 0
thermalBC(model, 'Edge', 3, 'Temperature', 0);

% Edge 4: Right boundary (z = 9, r = 2 → 4): u = 0.5 * (4/3) * (r^2 - 2^2)
thermalBC(model, 'Edge', 4, 'Temperature', @(region, state) ...
    0.5 * (4/3) * (region.y.^2 - 2^2));

% Edge 5: Upper boundary (z = 9 → 0, r = 4): u = 8
thermalBC(model, 'Edge', 5, 'Temperature', 8);

% Edge 6: Left boundary (z = 0, r = 4 → 0): u = 0.5 * r^2
thermalBC(model, 'Edge', 6, 'Temperature', @(region, state) ...
   0.5 * region.y.^2);

% Define thermal properties for Laplace's equation
thermalProperties(model, 'ThermalConductivity', 1);

% Generate the mesh
generateMesh(model);

% Solve the PDE
result = solve(model);

% Extract the solution
u = result.Temperature;

% Define the query point in the (z, r) plane
queryPoint = [3, 2];  % Replace zQuery and rQuery with the coordinates


% Interpolate the solution at the query point
uAtPoint = interpolateTemperature(result, queryPoint(1), queryPoint(2));

% Display the solution at the specified point
fprintf('The solution at (z = %.2f, r = %.2f) is %.4f\n', ...
        queryPoint(1), queryPoint(2), uAtPoint);