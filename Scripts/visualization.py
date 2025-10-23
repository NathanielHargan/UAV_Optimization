import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D

def drone_geometry(drone):
    ax = plt.figure().add_subplot(projection='3d')
    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-5, 5])
    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    for beam in drone_geo.beams:
        ax.plot([beam.start_node.location[0], beam.end_node.location[0]],
                [beam.start_node.location[1], beam.end_node.location[1]],
                [beam.start_node.location[2], beam.end_node.location[2]])


def drone_geometry_deformation(drone, m, q):
    # m: deformation multiplier
    # q: quality of image
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plt.title("UAV Slice Deformation")

    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-500, 500])

    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_zlabel("Z-axis (mm)")

    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    # Find the maximum deformation
    mag_max = 0
    for beam in drone_geo.beams:
        mag_0 = np.linalg.norm([beam.x_displacements_global[0],
                                beam.y_displacements_global[0],
                                beam.z_displacements_global[0]])
        mag_1 = np.linalg.norm([beam.x_displacements_global[1],
                                beam.y_displacements_global[1],
                                beam.z_displacements_global[1]])
        mag_max = max(mag_max, mag_0, mag_1)

    segments = []
    colors = []

    for beam in drone_geo.beams:
        deform_x = beam.x_displacements_global[0]
        deform_y = beam.y_displacements_global[0]
        deform_z = beam.z_displacements_global[0]
        position = beam.start_node.location
        for n in np.arange(0, 1, 1/q):
            deform_vector = beam.return_deformation_global(n + 1/q)
            position_n =  beam.start_node.location + (n+1/q) * (beam.end_node.location - beam.start_node.location)

            deform_x_n = deform_vector[0][0]
            deform_y_n = deform_vector[1][0]
            deform_z_n = deform_vector[2][0]
            deform_mag = np.sqrt(deform_x_n ** 2 + deform_y_n ** 2 + deform_z_n ** 2)/2 + np.sqrt(deform_x ** 2 + deform_y ** 2 + deform_z ** 2)/2

            col = mpl.cm.plasma(deform_mag/mag_max)
            segments.append([[position[0] + deform_x * m, position[1] + deform_y * m, position[2] + deform_z * m],
                            [position_n[0] + deform_x_n * m, position_n[1] + deform_y_n * m, position_n[2] + deform_z_n * m]])

            colors.append(col)
            position = position_n
            deform_x = deform_x_n
            deform_y = deform_y_n
            deform_z = deform_z_n

    # Create a Line3DCollection
    line_collection = Line3DCollection(segments, colors=colors, linewidth=2)
    ax.add_collection(line_collection)

    # Add colorbar
    sm = mpl.cm.ScalarMappable(cmap=mpl.cm.plasma, norm=mpl.colors.Normalize(vmin=0, vmax=mag_max))
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label("Deformation Magnitude (mm)")

    plt.show()


def drone_geometry_stress(drone, m, q):
    # m: deformation multiplier
    # q: quality of image
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plt.title("UAV Slice Bending Stress")

    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-500, 500])

    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_zlabel("Z-axis (mm)")

    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    # Find the maximum deformation
    stress_max = 0
    for beam in drone_geo.beams:
        max_0 = beam.bending_stresses_0
        max_1 = beam.bending_stresses_1
        stress_max = max(stress_max, max_0, max_1)


    stress_min = 0
    for beam in drone_geo.beams:
        min_0 = beam.bending_stresses_0
        min_1 = beam.bending_stresses_1
        stress_min = min(stress_min, min_0, min_1)

    segments = []
    colors = []

    for beam in drone_geo.beams:
        deform_x = beam.x_displacements_global[0]
        deform_y = beam.y_displacements_global[0]
        deform_z = beam.z_displacements_global[0]
        position = beam.start_node.location
        for n in np.arange(0, 1, 1/q):
            deform_vector = beam.return_deformation_global(n + 1/q)
            position_n =  beam.start_node.location + (n+1/q) * (beam.end_node.location - beam.start_node.location)

            deform_x_n = deform_vector[0][0]
            deform_y_n = deform_vector[1][0]
            deform_z_n = deform_vector[2][0]

            strain_bending, stress_bending, transverse_shear, transverse_shear_stress = beam.strain_stress_at_point_circle(n + 0.5/q)
            col = mpl.cm.plasma((stress_bending - stress_min)/(stress_max - stress_min))
            segments.append([[position[0] + deform_x * m, position[1] + deform_y * m, position[2] + deform_z * m],
                             [position_n[0] + deform_x_n * m, position_n[1] + deform_y_n * m, position_n[2] + deform_z_n * m]])

            colors.append(col)
            position = position_n
            deform_x = deform_x_n
            deform_y = deform_y_n
            deform_z = deform_z_n

    # Create a Line3DCollection
    line_collection = Line3DCollection(segments, colors=colors, linewidth=2)
    ax.add_collection(line_collection)

    # Add colorbar
    sm = mpl.cm.ScalarMappable(cmap=mpl.cm.plasma, norm=mpl.colors.Normalize(vmin=stress_min, vmax=stress_max))
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label("Bending Stress (MPa)")

    plt.show()


def drone_geometry_strain(drone, m, q):
    # m: deformation multiplier
    # q: quality of image
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plt.title("UAV Slice Bending Strain")

    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-500, 500])

    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_zlabel("Z-axis (mm)")

    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    # Find the maximum deformation
    strain_max = 0
    for beam in drone_geo.beams:
        max_0 = beam.bending_strains_0
        max_1 = beam.bending_strains_1
        strain_max = max(strain_max, max_0, max_1)

    strain_min = 0
    for beam in drone_geo.beams:
        min_0 = beam.bending_strains_0
        min_1 = beam.bending_strains_1
        strain_min = min(strain_max, min_0, min_1)


    segments = []
    colors = []

    for beam in drone_geo.beams:
        deform_x = beam.x_displacements_global[0]
        deform_y = beam.y_displacements_global[0]
        deform_z = beam.z_displacements_global[0]
        position = beam.start_node.location
        for n in np.arange(0, 1, 1/q):
            deform_vector = beam.return_deformation_global(n + 1/q)
            position_n = beam.start_node.location + (n+1/q) * (beam.end_node.location - beam.start_node.location)

            deform_x_n = deform_vector[0][0]
            deform_y_n = deform_vector[1][0]
            deform_z_n = deform_vector[2][0]

            strain_bending, stress_bending, transverse_shear, transverse_shear_stress = beam.strain_stress_at_point_circle(n + 0.5/q)
            col = mpl.cm.plasma((strain_bending - strain_min)/(strain_max - strain_min))
            segments.append([[position[0] + deform_x * m, position[1] + deform_y * m, position[2] + deform_z * m],
                             [position_n[0] + deform_x_n * m, position_n[1] + deform_y_n * m, position_n[2] + deform_z_n * m]])

            colors.append(col)
            position = position_n
            deform_x = deform_x_n
            deform_y = deform_y_n
            deform_z = deform_z_n

    # Create a Line3DCollection
    line_collection = Line3DCollection(segments, colors=colors, linewidth=2)
    ax.add_collection(line_collection)

    # Add colorbar
    sm = mpl.cm.ScalarMappable(cmap=mpl.cm.plasma, norm=mpl.colors.Normalize(vmin=strain_min, vmax=strain_max))
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label("Bending Strain (mm/mm)")

    plt.show()


def drone_geometry_strain_transverse_shear(drone, m, q):
    # m: deformation multiplier
    # q: quality of image
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plt.title("UAV Slice Transverse Shear Strain")

    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-500, 500])

    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_zlabel("Z-axis (mm)")

    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    # Find the maximum deformation
    transverse_shear_max = 0
    for beam in drone_geo.beams:
        max_0 = beam.strain_transverse_shear_0
        max_1 = beam.strain_transverse_shear_1
        transverse_shear_max = max(transverse_shear_max, max_0, max_1)


    transverse_shear_min = 0
    for beam in drone_geo.beams:
        min_0 = beam.strain_transverse_shear_0
        min_1 = beam.strain_transverse_shear_1
        transverse_shear_min = min(transverse_shear_min, min_0, min_1)

    segments = []
    colors = []

    for beam in drone_geo.beams:
        deform_x = beam.x_displacements_global[0]
        deform_y = beam.y_displacements_global[0]
        deform_z = beam.z_displacements_global[0]
        position = beam.start_node.location
        for n in np.arange(0, 1, 1/q):
            deform_vector = beam.return_deformation_global(n + 1/q)
            position_n =  beam.start_node.location + (n+1/q) * (beam.end_node.location - beam.start_node.location)

            deform_x_n = deform_vector[0][0]
            deform_y_n = deform_vector[1][0]
            deform_z_n = deform_vector[2][0]

            strain_bending, stress_bending, transverse_shear, transverse_shear_stress = beam.strain_stress_at_point_circle(n + 0.5/q)
            col = mpl.cm.plasma((transverse_shear - transverse_shear_min)/(transverse_shear_max - transverse_shear_min))
            segments.append([[position[0] + deform_x * m, position[1] + deform_y * m, position[2] + deform_z * m],
                             [position_n[0] + deform_x_n * m, position_n[1] + deform_y_n * m, position_n[2] + deform_z_n * m]])

            colors.append(col)
            position = position_n
            deform_x = deform_x_n
            deform_y = deform_y_n
            deform_z = deform_z_n

    # Create a Line3DCollection
    line_collection = Line3DCollection(segments, colors=colors, linewidth=2)
    ax.add_collection(line_collection)

    # Add colorbar
    sm = mpl.cm.ScalarMappable(cmap=mpl.cm.plasma, norm=mpl.colors.Normalize(vmin=transverse_shear_min, vmax=transverse_shear_max))
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label("Transverse Shear Strain (mm/mm)")

    plt.show()


def drone_geometry_stress_transverse_shear(drone, m, q):
    # m: deformation multiplier
    # q: quality of image
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plt.title("UAV Transverse Shear Stress")

    ax.set_xlim([0, 1000])
    ax.set_ylim([-500, 500])
    ax.set_zlim([-500, 500])

    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_zlabel("Z-axis (mm)")

    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)

    # Find the maximum deformation
    stress_transverse_shear_max = 0
    for beam in drone_geo.beams:
        max_0 = beam.stress_transverse_shear_0
        max_1 = beam.stress_transverse_shear_1
        stress_transverse_shear_max = max(stress_transverse_shear_max, max_0, max_1)


    stress_transverse_shear_min = 0
    for beam in drone_geo.beams:
        min_0 = beam.stress_transverse_shear_0
        min_1 = beam.stress_transverse_shear_1
        stress_transverse_shear_min = min(stress_transverse_shear_min, min_0, min_1)

    segments = []
    colors = []

    for beam in drone_geo.beams:
        deform_x = beam.x_displacements_global[0]
        deform_y = beam.y_displacements_global[0]
        deform_z = beam.z_displacements_global[0]
        position = beam.start_node.location
        for n in np.arange(0, 1, 1/q):
            deform_vector = beam.return_deformation_global(n + 1/q)
            position_n =  beam.start_node.location + (n+1/q) * (beam.end_node.location - beam.start_node.location)

            deform_x_n = deform_vector[0][0]
            deform_y_n = deform_vector[1][0]
            deform_z_n = deform_vector[2][0]

            strain_bending, stress_bending, transverse_shear, stress_transverse_shear = beam.strain_stress_at_point_circle(n + 0.5/q)
            col = mpl.cm.plasma((stress_transverse_shear - stress_transverse_shear_min) / (stress_transverse_shear_max - stress_transverse_shear_min))
            segments.append([[position[0] + deform_x * m, position[1] + deform_y * m, position[2] + deform_z * m],
                             [position_n[0] + deform_x_n * m, position_n[1] + deform_y_n * m, position_n[2] + deform_z_n * m]])

            colors.append(col)
            position = position_n
            deform_x = deform_x_n
            deform_y = deform_y_n
            deform_z = deform_z_n

    # Create a Line3DCollection
    line_collection = Line3DCollection(segments, colors=colors, linewidth=2)
    ax.add_collection(line_collection)

    # Add colorbar
    sm = mpl.cm.ScalarMappable(cmap=mpl.cm.plasma, norm=mpl.colors.Normalize(vmin=stress_transverse_shear_min, vmax=stress_transverse_shear_max))
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label("Transverse Shear Stress (MPa)")

    plt.show()


def design_space_mass(opt, design_variable_x, design_variable_y, sample_num, design_variables=None, x_var_bounds=None, y_var_bounds=None, active_constraints=None):
    # opt is opt object
    # constraints is an array of booleans
    if x_var_bounds is None:
        x_var_bounds = opt.boundaries[design_variable_x]
    if  y_var_bounds is None:
        y_var_bounds = opt.boundaries[design_variable_y]

    x_space = np.linspace(x_var_bounds[0], x_var_bounds[1], sample_num)
    y_space = np.linspace(y_var_bounds[0], y_var_bounds[1], sample_num)

    if design_variables is None:
        design_variables = opt.design_variables_initial_guess

    if active_constraints is None:
        active_constraints = opt.active_constraints_defaults

    cons_num = sum(1 for value in active_constraints.values() if value == True)

    constraints_map = np.zeros((sample_num, sample_num, cons_num))
    var_map = np.zeros((sample_num, sample_num))
    labels = []

    for x_i, x in enumerate(x_space):
        for y_i, y in enumerate(y_space):
            design_variables[design_variable_x] = x
            design_variables[design_variable_y] = y
            constraints_res, labels = opt.constraint_calculations(design_variables, active_constraints)
            constraints_map[x_i, y_i] = constraints_res
            var_map[x_i, y_i] = opt.mass(design_variables) * 1000
            
    plt.figure(figsize=(10, 10))

    img = plt.imshow(var_map, cmap='turbo', interpolation='nearest',
                     extent=[x_var_bounds[0], x_var_bounds[1], y_var_bounds[1], y_var_bounds[0]])
    cbar = plt.colorbar(img)

    # Set custom tick labels
    cbar.set_ticks(np.linspace(np.min(var_map), np.max(var_map), num=5))

    cbar.set_label('Mass (kg)')

    colors = ['r', 'g', 'b', 'y', 'm', 'c', 'k', 'orange', 'purple', 'brown', 'pink', 'gray']
    lvls = np.array([0.95, 1, 1.05])
    contours = []
    for i in range(cons_num):
        contours.append(plt.contour(x_space, y_space, constraints_map[:, :, i], levels=lvls, colors=colors[i]))

    for i in range(cons_num):
        plt.clabel(contours[i], lvls, inline=True, fontsize=8, fmt='%1.3f')
        plt.gca().set_aspect('auto', adjustable='box')

    legend_elements = []
    for i in range(cons_num):
        legend_elements.append(Line2D([0], [0], lw=2, label=labels[i], color=colors[i]))

    # Add the legend to the plot
    plt.legend(handles=legend_elements, loc='upper right')


def design_space_freq(opt, design_variable_x, design_variable_y, sample_num, design_variables=None, x_var_bounds=None,
                      y_var_bounds=None, active_constraints=None):
    # opt is opt object
    # constraints is an array of booleans
    if x_var_bounds is None:
        x_var_bounds = opt.boundaries[design_variable_x]
    if y_var_bounds is None:
        y_var_bounds = opt.boundaries[design_variable_y]

    x_space = np.linspace(x_var_bounds[0], x_var_bounds[1], sample_num)
    y_space = np.linspace(y_var_bounds[0], y_var_bounds[1], sample_num)

    if design_variables is None:
        design_variables = opt.design_variables_initial_guess

    if active_constraints is None:
        active_constraints = opt.active_constraints_defaults

    cons_num = sum(1 for value in active_constraints.values() if value == True)

    constraints_map = np.zeros((sample_num, sample_num, cons_num))
    var_map = np.zeros((sample_num, sample_num))
    labels = []

    for x_i, x in enumerate(x_space):
        for y_i, y in enumerate(y_space):
            design_variables[design_variable_x] = x
            design_variables[design_variable_y] = y
            constraints_res, labels = opt.constraint_calculations(design_variables, active_constraints, True)
            constraints_map[x_i, y_i] = constraints_res
            var_map[x_i, y_i] = opt.max_amp

    plt.figure(figsize=(10, 10))

    img = plt.imshow(var_map, cmap='turbo', interpolation='nearest',
                     extent=[x_var_bounds[0], x_var_bounds[1], y_var_bounds[1], y_var_bounds[0]],
                     norm=mpl.colors.LogNorm(vmin=var_map.min(), vmax=var_map.max()))

    cbar = plt.colorbar(img)

    # Set custom tick labels
    cbar.set_ticks(np.logspace(np.log10(var_map.min()), np.log10(var_map.max()), num=5))

    cbar.set_label('Amplitude (mm)')

    colors = ['y', 'm', 'c', 'k', 'orange', 'purple', 'brown', 'pink', 'gray']
    lvls = np.array([0.95, 1, 1.05])
    contours = []
    for i in range(cons_num):
        contours.append(plt.contour(x_space, y_space, constraints_map[:, :, i], levels=lvls, colors=colors[i]))

    for i in range(cons_num):
        plt.clabel(contours[i], lvls, inline=True, fontsize=8, fmt='%1.3f')
        plt.gca().set_aspect('auto', adjustable='box')

    legend_elements = []
    for i in range(cons_num):
        legend_elements.append(Line2D([0], [0], lw=2, label=labels[i], color=colors[i]))

    # Add the legend to the plot
    plt.legend(handles=legend_elements, loc='upper right')


def plot_drone_geometry(drone):
    ax = plt.figure(figsize=(15, 15)).add_subplot(projection='3d')
    ax.set_xlim([-1000, 1000])
    ax.set_ylim([-1000, 1000])
    ax.set_zlim([-5, 5])
    drone_geo = drone.beam_system
    node_loc = []
    for i in range(len(drone_geo.nodes)):
        node_loc.append(drone_geo.nodes[i].location)
    # ax.scatter(node_loc[:][0] + drone_geo.x_displacements,
    #           node_loc[:][1] + drone_geo.y_displacements,
    #           node_loc[:][2] + drone_geo.z_displacements)

    for beam in drone_geo.beams:
        ax.plot([beam.start_node.location[0], beam.end_node.location[0]],
                [beam.start_node.location[1], beam.end_node.location[1]],
                [beam.start_node.location[2], beam.end_node.location[2]], color='k')