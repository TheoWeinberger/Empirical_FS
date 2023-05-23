"""
Copyright 2023 T I Weinberger

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Code to empirically fit QO data to a Fermi surface defined by 
a geometric shape as opposed to band structure calculations

Note the code currently doesn't perform any symmetry operations/checking
and so this must be imposed by hand.

The code is currently set up for the calculation of QOs from a minimally corrugated
UTe2 fermi surface. This means that the symmetries used are for the case of Immm UTe2
and any other user must check the code carefully to remove undesired symmetry
operations.

Furthermore, current the code to extract extremal frequencies will only extract the 
maximum and minimum frequencies. This is because, in the case of the UTe2 cylinders, there
is no z-dependent warping of the cross-sectional area meaning that there is a maximum 
of two frequencies per sheet. This code must be adapted for more complex FS geometries
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.linalg import norm
import pyvista as pv
import triangle as tr
import pandas as pd
import seaborn as sns
import copy

#define some plotting parameters for rough
#plotting of output 
plt.style.use(['nature'])
cp = sns.color_palette("Set2")
pv.rcParams['transparent_background'] = True


def d_cylinder_1(p0, p1, Ra, Rb, theta_res):
    """
    Function to define the hole cylinder

    Args:
        p0: The first coordinate along the centerline of the cylinder
        p1The second coordinate along the centerline of the cylinder
        Ra: The radius of the cylinder in the a direction
        Rb: The radius of the cylinder in the b direction
        theta_res: number of angle intervals at which the surface is sampled

    Returns:
        X: X coordinates of the surface
        Y: Y coordinates of the surface
        Z: Z coordinates of the surface
    """

    #axis and radius
    p0 = np.array(p0)
    p1 = np.array(p1)

    #define a vector in direction of propogation
    v = np.array([0,0,0.01])

    #find magnitude of vector
    mag = norm(v)

    #normalise the vector
    v = v / mag
    mag = norm(p0-p1)

    #now define an abitrary vector that is normal to the direction
    #of propogation to define the full basis for the system
    v_perp = np.array([1, 0, 0])
    if (v == v_perp).all():
        v_perp = np.array([0, 1, 0])

    #make normal vector perpendicular to v
    n1 = np.cross(v, v_perp)

    #normalize n1
    n1 /= norm(n1)

    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.array([0])
    theta = np.linspace(0, 2*np.pi, theta_res)

    #use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)

    #generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + Rb * np.power(np.abs(np.sin(theta)),1/5)*np.sign(np.sin(theta)) * n1[i] + Ra * np.power(np.abs(np.cos(theta)),1/5)*np.sign(np.cos(theta)) *n2[i]  for i in [0, 1, 2]]
    X = X*np.power(1-0.6*np.cos(4*theta), 0.05)
    Y = Y*np.power(1-0.6*np.cos(4*theta), 0.05)

    return X,Y,Z


def d_cylinder_2(p0, p1, Ra, Rb, theta_res):
    """
    Function to define the hole cylinder

    Args:
        p0: The first coordinate along the centerline of the cylinder
        p1The second coordinate along the centerline of the cylinder
        Ra: The radius of the cylinder in the a direction
        Rb: The radius of the cylinder in the b direction
        theta_res: number of angle intervals at which the surface is sampled

    Returns:
        X: X coordinates of the surface
        Y: Y coordinates of the surface
        Z: Z coordinates of the surface
    """

    #axis and radius
    p0 = np.array(p0)
    p1 = np.array(p1)

    #define a vector in direction of propogation
    v = np.array([0,0,0.01])

    #find magnitude of vector
    mag = norm(v)

    #normalise the vector
    v = v / mag
    mag = norm(p0-p1)

    #now define an abitrary vector that is normal to the direction
    #of propogation to define the full basis for the system
    v_perp = np.array([1, 0, 0])
    if (v == v_perp).all():
        v_perp = np.array([0, 1, 0])

    #make normal vector perpendicular to v
    n1 = np.cross(v, v_perp)

    #normalize n1
    n1 /= norm(n1)

    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.array([0])
    theta = np.linspace(0, 2*np.pi, theta_res)

    #use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)

    #generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + Rb * np.power(np.abs(np.sin(theta)),1/5)*np.sign(np.sin(theta)) * n1[i] + Ra * np.power(np.abs(np.cos(theta)),1/5)*np.sign(np.cos(theta)) *n2[i]  for i in [0, 1, 2]]
    Y = Y*(1+np.power(1+0.8*np.cos(p0[2])*np.cos((theta)), 0.5))
    Z = Z

    return X,Y,Z


def create_surface(surface_func, r, t, Ra, Rb, theta_res):
    """
    Function to create a surface given a function, for generalisation 
    purposed this should be redefined to accept **kwargs

    Args:
        surface_func: the function which defines the surface
        r: the degree of warping
        t: the number of iterations of the surface in the z direction
        Ra: the extent of the surface in the a/x-direction
        Rb: the extent of the surface in the b/y-direction
        theta_res: number of angle intervals at which the surface is sampled

    Returns:
        X: X coordinates of the surface
        Y: Y coordinates of the surface
        Z: Z coordinates of the surface
    """


    u = np.linspace(-t*np.pi, t*np.pi, res)
    x = r*np.cos(u)
    y = np.zeros_like(u)
    z = u


    coord1 = [x[0], y[0], z[0]]
    coord2 = [x[1], y[1], z[1]]

    X,Y,Z = surface_func(coord1, coord2, Ra, Rb, theta_res)

    for i in range(1, len(u) - 1):
        coord1 = [x[i], y[i], z[i]]
        coord2 = [x[i+1], y[i+1], z[i+1]]

        X_temp,Y_temp,Z_temp = surface_func(coord1, coord2, Ra, Rb, theta_res)

        X = np.concatenate((X, X_temp))
        Y = np.concatenate((Y, Y_temp))
        Z = np.concatenate((Z, Z_temp))

    return X,Y,Z


def triangulate_slice(slc):
    """
    Code to triangulate the slice taking through the cylinders
    when calculating extremal frequencies. This step is required
    to turn the outline of the slice through the fermi surface
    into a 2D shape from which the area can be caluclated

    Args:

        slc: the slice through the fermi surface

    Returns:

        pd: the triangulated slice from which the area can be determined
    """

    seg = slc.lines.reshape(-1, 3)[:, 1:]
    A = dict(vertices=np.array(slc.points[:, :2], dtype=np.float64), segments=seg)
    B = tr.triangulate(A, 'pci')

    n_faces = B['triangles'].shape[0]
    triangles = np.hstack((np.ones((n_faces, 1), np.int32)*3, B['triangles']))

    # back to 3D
    pts = np.empty((B['vertices'].shape[0], 3))
    pts[:, :2] = B['vertices']
    pts[:, -1] = slc.points[0, 2]
    pd = pv.PolyData(pts, triangles, n_faces)

    return pd


def get_brillouin_zone_3d(cell):
    """
    Uses the k-space vectors and voronoi analysis to define 
    the BZ of the system

    Args:
        cell: a 3x3 matrix defining the basis vectors in 
        reciprocal space

    Returns:
        vor.vertices[bz_vertices]: vertices of BZ
        bz_ridges: edges of the BZ
        bz_facets: BZ facets

    """


    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):

        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))
    
    return vor.vertices[bz_vertices], bz_ridges, bz_facets


def calculate_extremal_freqs(direction, surfaces, i, norm, plane_res, rotation_points, theta_res):
    """
    Calculate the maximum and minimum frequencies associated with a fermi surface sheet
    at each field angle. In the case of UTe2 here, this is valid since there is no
    z-direction modulation in the FS area and hence there are a maximum of two different
    frequencies at each angle

    Args:
        direction: the direciton through which the field is rotated, either c-a or c-b
        surfaces: list of the surfaces
        i:the index of the surface
        norm: the normalisation of the fermi surface used to calculate real frequencies
        plane_res: the number of planes along the z access for sampling the area
        rotation_points: the number of angles to be sampled as the field rotates through 90 degrees
        theta_res number of angle intervals at which the surface is sampled
    """
        
    theta = np.linspace(0,np.pi/2,rotation_points)

    origins = [(0,0,2*np.pi*(i-plane_res/2)/(plane_res/2)) for i in range(plane_res)]

    grid = pv.StructuredGrid()
    grid.points = np.array([surfaces[0].ravel(), surfaces[1].ravel(), surfaces[2].ravel()]).T
    grid.dimensions = [1, theta_res[i], res - 1]
    grid = grid.extract_surface()

    max_area = []
    min_area = []
    print("Analysing Surface: " + str(i))
    print("Direction: " + direction)
    for theta_val in theta:
        print("Theta: " + str(theta_val*360/(2.0*np.pi)))
        areas = []
        for origin in origins:
            slice_z = 0
            if direction == 'c-a':
                slice_z = grid.slice(normal=(np.sin(theta_val),0,np.cos(theta_val)), origin=origin, generate_triangles=True)
                slice_z.rotate_y(-theta_val*360/(2.0*np.pi), inplace=True)
            elif direction == 'c-b':
                slice_z = grid.slice(normal=(0,np.sin(theta_val),np.cos(theta_val)), origin=origin, generate_triangles=True)
                slice_z.rotate_x(theta_val*360/(2.0*np.pi), inplace=True)
            closed_slice = triangulate_slice(slice_z)
            areas += [np.sum(closed_slice.compute_cell_sizes()["Area"])]

        max_area += [max(areas)]
        min_area += [min(areas)]

    renorm = (37409.6788)/((norm/(2*np.pi))**2)
    renorm = renorm/1000
    max_area = [area*renorm for area in max_area]
    min_area = [area*renorm for area in min_area]

    df = pd.DataFrame (theta*360/(2.0*np.pi), columns = ['angle'])
    df["max_freq"] = max_area
    df["min_freq"] = min_area

    df.to_csv("UTe2_band" + str(i) + "." + direction + ".out", sep=',', header=True, index=False)


## Input Parameters for hole sheet
n_sheets = 2
stretch = 1.05
res = 2000   # plot resolution; higher -> finer resolution
theta_res= 100
Ra_h = 2.02/stretch   # squircle radius in the a direction
Rb_h = 2.02*(stretch) # squircle radius in the b direction
r_h = 0.52      # degree of warping
t_h = 20      # number of BZ in z direction

## Input Parameters for electron sheet
Ra_e = 1.295 # squircle radius in the a direction (for the electron sheet this value is half the radius)
Rb_e = 1.55 # squircle radius in the b direction
r_e = -0.15  # degree of warping

##input parameters for frequency analysis:
plane_res = 60 # number of planes along the B field direction for intersection analysis
rotation_points = 45 #the number of angles for frequency calculations


# combine all parameters into lists for iterating over multiple sheets
Ra = [Ra_h, Ra_e]
Rb = [Rb_h, Rb_e]
r = [r_h, r_e]
t = [t_h, t_h]
theta_res = [theta_res, theta_res]
func = [d_cylinder_1, d_cylinder_2]
lattice_parameters = [7.791341, 11.500873, 26.100897] #in bohr radii


"""
Below here is where the real calculation occurs
"""
surfaces = []

for i in range(n_sheets):

    # Create the fermi surface for the hole sheet
    #define the centerline of the cylinder
    X, Y, Z = create_surface(func[i], r[i], t[i], Ra[i], Rb[i], theta_res[i])
    surfaces += [[X,Y,Z]]

#because of the way the electron sheet is parameterised it needs to be rotated 90 degrees
#########################################################################################
#Generally speaking a surface should be defined to not require a rotation therefore any other
#users should comment this out
#########################################################################################
surface_temp = copy.deepcopy(surfaces)
surfaces[1][0] = surface_temp[1][1]
surfaces[1][1] = surface_temp[1][0]

k = [1/lattice_parameter for lattice_parameter in lattice_parameters]

norm = np.pi/(k[2])

kx = k[0]*norm
ky = k[1]*norm
kz = k[2]*norm

#define reciprocal space lattice vectors for UTe2
#################################################
#again this is specific to UTe2/body centred systems
#################################################

cell = np.array([[0, ky, kz],
                 [kx, 0, kz],
                 [kx, ky, 0]])

plotter = pv.Plotter(lighting="three lights")

#generate BZ from voronoi analysis
v, e, f = get_brillouin_zone_3d(cell)

#triangulate BZ surface
bz_surf = pv.PolyData(v)
bz_surf = bz_surf.delaunay_3d()
bz_surf = bz_surf.extract_surface()
edges = bz_surf.extract_all_edges()

for i in range(n_sheets):

    if i == 1:
        
        grid1 = pv.StructuredGrid()
        grid1.points = np.array([surfaces[i][0].ravel(), (surfaces[i][1] + ky/2).ravel(), surfaces[i][2].ravel()]).T
        grid1.dimensions = [1, theta_res[i], res - 1]
        grid1 = grid1.extract_surface()

        grid2 = pv.StructuredGrid()
        grid2.points = np.array([-surfaces[i][0].ravel(), (-surfaces[i][1] - ky/2).ravel(), surfaces[i][2].ravel()]).T
        grid2.dimensions = [1, theta_res[i], res - 1]
        grid2 = grid2.extract_surface()
    
        plotter.add_mesh(grid1.clip_surface(bz_surf, invert=True), lighting=True, color="#73D2DE", opacity=1.0)
        plotter.add_mesh(grid2.clip_surface(bz_surf, invert=True), lighting=True, color="#73D2DE", opacity=1.0)

    if i == 0:

        grid1 = pv.StructuredGrid()
        grid1.points = np.array([(surfaces[i][0] + kx/2).ravel(), -surfaces[i][1].ravel(), surfaces[i][2].ravel()]).T
        grid1.dimensions = [1, theta_res[i], res - 1]
        grid1 = grid1.extract_surface()

        grid2 = pv.StructuredGrid()
        grid2.points = np.array([(-surfaces[i][0] - kx/2).ravel(), surfaces[i][1].ravel(), surfaces[i][2].ravel()]).T
        grid2.dimensions = [1, theta_res[i], res - 1]
        grid2 = grid2.extract_surface()

        #for some reason clipping doesn't work properly on this surface
        #this part is a patch
        grid1.compute_implicit_distance(bz_surf, inplace=True)
        grid2.compute_implicit_distance(bz_surf, inplace=True)

        plane1 = np.array([kx, 0, kz])
        plane2 = np.array([kx, 0, -kz])
        plane3 = np.array([-kx, 0, -kz])
        plane4 = np.array([-kx, 0, kz])

        grid1 = grid1.threshold(0.0, scalars="implicit_distance", invert=True)
        grid1 = grid1.clip(normal= plane1, origin = plane1/2)
        grid1 = grid1.clip(normal= plane2, origin = plane2/2)

        grid2 = grid2.clip(normal= plane3, origin = plane3/2)
        grid2 = grid2.clip(normal= plane4, origin = plane4/2)

        plotter.add_mesh(grid1, lighting=True, color="#f5b105", opacity=1.0)
        plotter.add_mesh(grid2.threshold(0.0, scalars="implicit_distance", invert=True), lighting=True, color="#f5b105", opacity=1.0)

plotter.set_background('white')
plotter.camera_position = 'yz'
plotter.camera.azimuth = 65
plotter.camera.elevation = 30
plotter.enable_parallel_projection()

#plot BZ
for xx in e:
    line = pv.MultipleLines(points = np.array([xx[:, 0], xx[:, 1], xx[:, 2]]).T)
    plotter.add_mesh(line, color = "black") 

plotter.show()

#Now calculate the extremal frequencies associated with each sheet
for i in range(n_sheets):

    calculate_extremal_freqs('c-a', surfaces[i], i, norm, plane_res, rotation_points, theta_res)
    calculate_extremal_freqs('c-b', surfaces[i], i, norm, plane_res, rotation_points, theta_res)