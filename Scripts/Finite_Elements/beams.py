#%%
import numpy as np
import math

from Scripts.Finite_Elements.Beam_FEA_System.beam_system import BeamSystem
from Scripts.Finite_Elements.Beam_FEA_System.beam_type import BeamType

import logging
from numpy.testing import assert_almost_equal


def cantilever_rectangle():
    # Relevant
    beam_system = BeamSystem("Cantilever-end load")
    arm_beam = BeamType("rectangle", [20, 6], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_force(np.array([0, 0, -4]), "point_1")
    beam_system.add_boundary_condition(0,"x","point_1")
    beam_system.add_boundary_condition(0,"y","point_1")
    beam_system.add_boundary_condition(0,"theta x","point_1")
    beam_system.add_boundary_condition(0,"theta z","point_1")

    beam_system.solve_FEA()
    beam_system.solve_static()
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started Cantilever—end load")
    print("calculated z displacement: " + str(y) + " | FEA z displacement: " + str(beam_system.z_displacements[1]))
    print("calculated y displacement: " + str(0) + " | FEA y displacement: " + str(beam_system.y_displacements[1]))
    print("calculated moment: " + str(-m) + " | FEA moment reaction: " + str(beam_system.y_moments[0]))
    print("calculated reaction: " + str(F) + " | FEA reaction z: " + str(beam_system.z_forces[0]))
    assert_almost_equal(y, beam_system.z_displacements[1], 3)
    assert_almost_equal(0, beam_system.y_displacements[1])
    assert_almost_equal(-dydx, beam_system.y_angles[1])
    assert_almost_equal(-m, beam_system.y_moments[0])
    assert_almost_equal(F, beam_system.z_forces[0])

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))


def center_multidim():
    beam_system = BeamSystem("Multidimensional Simple Support Center Load")
    arm_beam = BeamType("annulus", [25, 20], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([-1000, -1000, -1000]), "point_0")
    beam_system.add_node(np.array([0, 0, 0]), "point_1")
    beam_system.add_node(np.array([1000, 1000, 1000]), "point_2")
    beam_system.add_beam(arm_beam, "point_0", "point_1", np.array([0, -1, 1]), "beam_1")
    beam_system.add_beam(arm_beam, "point_1", "point_2", np.array([0, -1, 1]), "beam_2")
    beam_system.add_boundary_condition(0, "x", "point_0")
    beam_system.add_boundary_condition(0, "y", "point_0")
    beam_system.add_boundary_condition(0, "z", "point_0")
    beam_system.add_boundary_condition(0, "x", "point_2")
    beam_system.add_boundary_condition(0, "y", "point_2")
    beam_system.add_boundary_condition(0, "z", "point_2")
    beam_system.add_force(np.array([0, 1, -1]), "point_1")
    beam_system.solve_FEA()

    L = math.sqrt(3 * 2000 ** 2)
    I = (math.pi/64) * (50**4 - 40**4)
    E = 71700
    F = math.sqrt(2)
    # kuu
    # Z1  Z2
    # Z2  Z3
    y = -(F * (L**3)) / (48 * E * I)
    r = F/2
    dydx = -(F * L ** 2) / (16 * E * I)
    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info('Started Multidimensional Simple Support Center Load Test')
    logging.debug("FEA x: " + str(beam_system.x_displacements[1]))
    logging.debug("FEA y: " + str(beam_system.y_displacements[1]))
    logging.debug("FEA z: " + str(beam_system.z_displacements[1]))
    logging.debug("FEA x theta: " + str(beam_system.x_angles[1]))
    logging.debug("FEA y theta: " + str(beam_system.y_angles[1]))
    logging.debug("FEA z theta: " + str(beam_system.z_angles[1]))
    logging.debug("FEA x reaction: " + str(beam_system.x_forces[0]))
    logging.debug("FEA y reaction: " + str(beam_system.y_forces[0]))
    logging.debug("FEA z reaction: " + str(beam_system.z_forces[0]))

    total_fea_disp = math.sqrt(beam_system.x_displacements[1] ** 2 + beam_system.y_displacements[1] ** 2 + beam_system.z_displacements[1] ** 2)
    total_fea_reaction = math.sqrt(beam_system.x_forces[0] ** 2 + beam_system.y_forces[0] ** 2 + beam_system.z_forces[0] ** 2)

    logging.debug("Calculated total displacement: " + str(-y) + " | Total displacement FEA: " + str(total_fea_disp))
    logging.debug("Calculated total reaction: " + str(r) + " | Total reaction FEA: " + str(total_fea_reaction))
    assert_almost_equal(-y, total_fea_disp, 3)
    assert_almost_equal(r, total_fea_reaction)
    logging.info('Finished')


def incline_boundary_conditions():
    beam_system = BeamSystem("Incline Boundary Conditions")
    beam_12 = BeamType("circle", [27.63953195], "Aluminum7075-T6", "beam_12")
    beam_3 = BeamType("circle", [32.869128059], "Aluminum7075-T6", "beam_3")

    beam_12.material_properties['elastic modulus'] = 210000
    beam_3.material_properties['elastic modulus'] = 210000
    beam_12.cross_section_properties['transverse shear deflection constant y'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant z'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant y'] = 0
    beam_12.cross_section_properties['transverse shear deflection constant z'] = 0
    print(beam_12.cross_section_properties)

    beam_system.add_node(np.array([0, 0, 0]), "point_1")
    beam_system.add_node(np.array([0, 1000, 0]), "point_2")
    beam_system.add_node(np.array([1000, 1000, 0]), "point_3")
    beam_system.add_beam(beam_12, "point_1", "point_2", np.array([0, 0, 1]), "beam_1")
    beam_system.add_beam(beam_12, "point_2", "point_3", np.array([0, 0, 1]), "beam_2")
    beam_system.add_beam(beam_3, "point_1", "point_3", np.array([0, 0, 1]), "beam_3")

    beam_system.add_boundary_condition(0, "x", "point_1")
    beam_system.add_boundary_condition(0, "y", "point_1")
    beam_system.add_boundary_condition(0, "z", "point_1")
    beam_system.add_boundary_condition(0, "theta x", "point_1")
    beam_system.add_boundary_condition(0, "theta y", "point_1")
    beam_system.add_boundary_condition(0, "theta z", "point_1")

    beam_system.add_boundary_condition(0, "y", "point_2")
    beam_system.add_boundary_condition(0, "z", "point_2")
    beam_system.add_boundary_condition(0, "theta y", "point_2")
    beam_system.add_boundary_condition(0, "theta z", "point_2")


    beam_system.add_boundary_condition(0, "y", "point_3", np.array([0.707, 0.707, 0]), np.array([-0.707, 0.707, 0]))
    beam_system.add_boundary_condition(0, "z", "point_3")
    beam_system.add_boundary_condition(0, "theta x", "point_3")
    beam_system.add_boundary_condition(0, "theta y", "point_3")
    beam_system.add_boundary_condition(0, "theta z", "point_3")
    beam_system.add_force(np.array([1000000, 0, 0]), "point_2")
    beam_system.solve_FEA()
    np.set_printoptions(linewidth=400)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info('Started Incline Boundary Conditions')
    logging.debug("\nPRE-TRANSFORM 3:\n" + str(beam_system.global_stiffness_matrix))
    logging.debug("\nTRANSFORM 3:\n" + str(beam_system.local_bc_transform))
    logging.debug("\nPOST-TRANSFORM 3:\n" + str(beam_system.applied_stiffness_matrix))
    logging.debug("FEA 2y: " + str(beam_system.x_displacements[1]) + ' | Expected 11.91')
    # assert_almost_equal(beam_system.x_displacements[1], 11.91, 3)
    logging.debug("FEA 1x: " + str(beam_system.x_forces[0]) + ' | Expected -500000')
    logging.debug("FEA 1x: " + str(beam_system.x_forces[0]) + ' | Expected -500000')
    # assert_almost_equal(beam_system.x_forces[0], -500000, -4)
    logging.debug("FEA 1y: " + str(beam_system.y_forces[0]) + ' | Expected -500000')
    logging.debug("FEA 1y: " + str(beam_system.y_forces[0]) + ' | Expected -500000')
    # assert_almost_equal(beam_system.y_forces[0], -500000, -4)
    logging.debug("FEA 2y: " + str(beam_system.y_forces[1]) + ' | Expected 0')
    # assert_almost_equal(beam_system.y_forces[1], 0, -4)
    logging.info('Finished')


def reversed_cantilever_rectangle():
    # Relevant
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    beam_system = BeamSystem("Reverse Cantilever-end load")
    arm_beam = BeamType("rectangle", [20, 6], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_boundary_condition( -0.5816413736845836,"z", "point_1")


    beam_system.solve_FEA()

    print("\nFORCE MOMENT VECTOR: \n" + str(beam_system.force_moment_vector))

    print("KUU: " + str(beam_system.kuu))
    print("KPP: " + str(beam_system.kpp))

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started Reverse Cantilever—end load")
    print("calculated z displacement: " + str(y) + " | FEA z displacement: " + str(beam_system.z_displacements[1]))
    print("calculated reaction: " + str(F) + " | FEA reaction z: " + str(beam_system.z_forces[1]))
    print("stress: ", beam_system.beams[0].stresses_bending_z)
    print("expected stress: ", m/I)

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))


def cantilever_annulus():
    beam_system = BeamSystem("annulus Cantilever-end load")
    arm_beam = BeamType("hexagon", [20], "Aluminum7075-T6", "arm_beam")
    beam_system.add_node(np.array([0,0,0]),"origin")
    beam_system.add_node(np.array([500,0,0]),"point_1")
    beam_system.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system.add_boundary_condition(0,"x","origin")
    beam_system.add_boundary_condition(0,"y","origin")
    beam_system.add_boundary_condition(0,"z","origin")
    beam_system.add_boundary_condition(0,"theta x","origin")
    beam_system.add_boundary_condition(0,"theta y","origin")
    beam_system.add_boundary_condition(0,"theta z","origin")
    beam_system.add_force(np.array([0, 0, -4]), "point_1")

    beam_system.solve_FEA()

    beam_system_2 = BeamSystem("annulus Cantilever-end load")
    beam_system_2.add_node(np.array([0,0,0]),"origin")
    beam_system_2.add_node(np.array([500,0,0]),"point_1")
    beam_system_2.add_beam(arm_beam,"origin","point_1",np.array([0,0,1]),"beam")
    beam_system_2.add_boundary_condition(0,"x","origin")
    beam_system_2.add_boundary_condition(0,"y","origin")
    beam_system_2.add_boundary_condition(0,"z","origin")
    beam_system_2.add_boundary_condition(0,"theta x","origin")
    beam_system_2.add_boundary_condition(0,"theta y","origin")
    beam_system_2.add_boundary_condition(0,"theta z","origin")
    beam_system_2.add_boundary_condition(-0.00917048964187906,"z","point_1")

    beam_system_2.solve_FEA()
    L = 500
    I = (6 * 20 ** 3) / 12
    E = 71700
    F = 4

    y = - (F * L ** 3) / (3*E*I)
    m = F * L
    dydx = (F*L**2)/(2*E*I) - (F*L**2)/(E*I)

    logging.basicConfig(filename='beams.log', level=logging.DEBUG)
    logging.info("Started annulus Cantilever—end load")
    print(" FEA z displacement: " + str(beam_system.z_displacements[1]))
    print(" FEA z force: " + str(beam_system_2.z_forces[1]))
    print("stress: ", beam_system.beams[0].stresses_bending_z)
    print("expected stress: ", m/beam_system.beams[0].beam_type.cross_section_properties['second moment of area y'])


    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape z disp " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][0]))

    for i in np.arange(0, 1.1, 0.1):
        logging.debug("Shape y angle " + str(round(i,2)) + ": " + str(beam_system.beams[0].return_shape_functions(i)[0][1]))


def beam_30_deg():
    beam_system = BeamSystem("30 Deg Axial Tension")
    arm_beam = BeamType("circle", [20], "Aluminum7075-T6", "arm_beam")

    # 2 length 200 elements
    beam_system.add_node(np.array([0, 0, 0]),"origin")
    beam_system.add_node(np.array([math.sqrt(3)*100, 100, 0]),"point_1")
    beam_system.add_node(np.array([500, 200, 0]),"point_2")

    beam_system.add_beam(arm_beam, "origin", "point_1", np.array([0, 0, 1]), "axial_tension_beam")
    beam_system.add_beam(arm_beam, "point_1", "point_2", np.array([0, 0, 1]), "pulling_beam")

    beam_system.add_boundary_condition(0, "x", "origin")
    beam_system.add_boundary_condition(0, "y", "origin")
    beam_system.add_boundary_condition(0, "z", "origin")
    beam_system.add_boundary_condition(0, "theta x", "origin")
    beam_system.add_boundary_condition(0, "theta y", "origin")
    beam_system.add_boundary_condition(0, "theta z", "origin")

    beam_system.add_boundary_condition(0, "y", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "z", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta x", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta y", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))
    beam_system.add_boundary_condition(0, "theta z", "point_1", np.array([math.sqrt(3), 1, 0]), np.array([0,0,1]))

    beam_system.add_force(np.array([0, 1, 0]), "point_2")

    beam_system.solve_FEA()
    beam_system.solve_static()

    N1_disp = beam_system.nodes[beam_system.select_node("point_1")].result_displacement
    print("Displacement 1: ", N1_disp)
    print("Displacement 2: ", beam_system.nodes[beam_system.select_node("point_2")].result_displacement)

    print("Displacement 1 theta: ", math.degrees(math.atan2(N1_disp[0], N1_disp[1])))


def timoshenko_simply_supp_point():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "y", node_count-1)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, 0]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920

    for beam in beam_system.beams:
        print("Max Stress 0", max(beam.bending_stresses_0))
        print("Max Stress 1", max(beam.bending_stresses_1))


def timoshenko_simply_supp():
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")

    node_count = 21
    L = 1000
    Q = 0.0222411 # Load
    E = arm_beam.material_properties['elastic modulus']
    I = arm_beam.cross_section_properties['second moment of area z']
    EI = E*I
    K = arm_beam.cross_section_properties['transverse shear deflection constant z']
    G = arm_beam.material_properties["shear modulus"]
    A = arm_beam.cross_section_properties['area']

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "x", 0)
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "x", node_count-1)
    beam_system.add_boundary_condition(0, "y", node_count-1)

    for i in range(node_count):
        beam_system.add_boundary_condition(0, "z symm", i)
    '''
    beam_system.add_force(np.array([0, -Q, 0]), 10)
    beam_system.solve_FEA()

    print("Midpoint Disp:", beam_system.nodes[10].result_displacement)
    print("Midpoint Rotation:", beam_system.nodes[10].result_displacement)
    # Abaqus Results
    #
    '''

    # Add distributed load
    q = L*Q/(node_count)  # load on each node
    for i in range(1, node_count-1):  # ingnoring phantom loads
        beam_system.add_force(np.array([0, -q, 0]), i)
        print("")

    beam_system.add_moment(np.array([0, 0, (Q * (L/(node_count-1)) ** 2)/12]), 0)
    beam_system.add_moment(np.array([0, 0, -(Q * (L/(node_count-1)) ** 2)/12]), node_count-1)

    beam_system.solve_FEA()
    xi_midpoint = int(node_count/2)
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[xi_midpoint].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    x = beam_system.nodes[xi_midpoint].location[0]
    print("x:", x)
    w_bending = (Q/EI) * ((x**4 / 24) - ((L**2 * x**2)/16) + (5 * L ** 4)/384)
    w_shear = (Q/(2 * K * G * A)) * (((L ** 2) / 4) - (x ** 2))
    w = w_bending + w_shear
    print("Expected Midpoint Disp from bending: ", w_bending)
    print("Expected Midpoint Disp from shear: ", w_shear)
    print("Expected Midpoint Disp from total: ", w)


    x = beam_system.nodes[xi_75].location[0]
    print("75:",x)
    w_bending = (Q/EI) * ((x**4 / 24) - ((L**2 * x**2)/16) + (5 * L ** 4)/384)
    w_shear = (Q/(2 * K * G * A)) * (((L ** 2) / 4) - (x ** 2))
    w = w_bending + w_shear
    print("Expected 75 Disp from bending: ", w_bending)
    print("Expected 75 Disp from shear: ", w_shear)
    print("Expected 75 Disp from total: ", w)


def timoshenko_simply_supp_point_shear():
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    L = 50
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x, 0, 0]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 1, 0]))

    # Create BCS
    beam_system.add_boundary_condition(0, "x", 0)
    beam_system.add_boundary_condition(0, "x", node_count-1)
    beam_system.add_boundary_condition(0, "y", 0)
    beam_system.add_boundary_condition(0, "y", node_count-1)

    for i in range(node_count):
        beam_system.add_boundary_condition(0, "z symm", i)

    print("~Point Load Shear Test~")
    beam_system.add_force(np.array([0, -Q, 0]), 10)

    beam_system.solve_FEA()

    xi_midpoint = int(node_count/2)
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[10].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -3.997e-05
    # ABAQUS 23 midpoint: -2.232e-05


def timoshenko_stress_circle():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("circle", [25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x+0.0001, 0, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, -Q]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920
    for beam in beam_system.beams:
        print("Max Stress 0", np.round(beam.bending_stresses_0,5))
        print("Max Stress 1", np.round(beam.bending_stresses_1,5))

    print(beam_system.y_angles)
    # print(beam_system.du)

    # print(beam_system.kuu)
    print(beam_system.z_angles)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.du)


def timoshenko_stress():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 21
    midpoint_index = int((node_count-1)/2)
    L = 1000
    Q = 20 # Load

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)


    print("~Point Load Slender~")
    beam_system.add_force(np.array([0, -Q, -Q]), midpoint_index)

    beam_system.solve_FEA()
    xi_midpoint = midpoint_index
    xi_75 = int(node_count*(3/4))
    print("Midpoint Disp:", beam_system.nodes[midpoint_index].result_displacement)
    print("0 rotation:", beam_system.nodes[0].result_angular_displacement)
    print("0.75 Disp:", beam_system.nodes[xi_75].result_displacement)

    # ABAQUS 22 midpoint: -0.178522mm
    # ABAQUS 23 midpoint: -0.178871mm
    # Abaqus max stress: 1.920
    for beam in beam_system.beams:
        print("Max Stress 0", np.round(beam.bending_stresses_0,5))
        print("Max Stress 1", np.round(beam.bending_stresses_1,5))

    # print(beam_system.y_angles)
    # print(beam_system.du)

    # print(beam_system.kuu)
    # print(beam_system.z_angles)
    print(np.shape(beam_system.global_stiffness_matrix))
    print(np.shape(beam_system.global_mass_matrix))

    beam_system.solve_natural_frequencies()
    print("nf:", beam_system.natural_frequencies)
    print("ev:", beam_system.natural_frequencies_eigenvals)
    print("ev shape:", np.shape(beam_system.natural_frequencies_eigenvals))
    print("k shape:", np.shape(beam_system.global_stiffness_matrix))

    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.beams[10].transformation_matrix)
    # print(beam_system.du)


def timoshenko_freq():
    #
    beam_system = BeamSystem("timoshenko_simply_supp")
    arm_beam = BeamType("rectangle", [25, 25], "Aluminum7075-T6", "arm_beam")
    node_count = 50
    midpoint_index = int((node_count-1)/2)
    L = 1000

    # Create nodes
    for i in range(node_count):
        x = L * (i/(node_count-1)) - (L/2)
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    # Create elements
    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    # Create BCS
    beam_system.add_boundary_condition(0, "pinned", 0)
    beam_system.add_boundary_condition(0, "pinned", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", node_count-1)
    beam_system.add_boundary_condition(0, "theta x", 0)
    beam_system.solve_FEA()

    np.set_printoptions(linewidth=4000)
    # print(np.shape(beam_system.kuu))
    # print(np.shape(beam_system.muu))

    # print("kuu,", np.round(beam_system.kuu,7))
    print("muu,", np.round(beam_system.muu,7))
    # print("m,", np.round(beam_system.global_mass_matrix,7))
    # print("m_e,", np.round(beam_system.beams[0].global_mass_matrix,7))

    beam_system.solve_natural_frequencies()
    # print("lamda:", beam_system.natural_frequencies_eigenvals)
    # print("omega:", beam_system.natural_frequencies)
    print("Hz:", np.round(beam_system.natural_frequencies/(2*math.pi)))
    # print("ev shape:", np.shape(beam_system.natural_frequencies_eigenvals))
    # print("k shape:", np.shape(beam_system.global_stiffness_matrix))

    # m = 8.0605e-4 * np.array([[2,0],[0,1]])
    # k = 2.6507e6 * np.array([[2,-1],[-1,1]])
    # print("m:",m)
    # print("k:",k)
    # eigenvals = scipy.linalg.eigvals(k, m)
    # print("e,", eigenvals)
    # print("w,", np.sqrt(eigenvals))
    # print("hz,", np.sqrt(eigenvals)/(2*math.pi))


def transient_displacement():
    timesteps = 100
    node_count = 3
    L = 100

    beam_system = BeamSystem("canilever beam")
    arm_beam = BeamType("annulus", [14, 12], "Carbon Fiber", "arm_beam")

    for i in range(node_count):
        x = L * (i/(node_count-1))
        beam_system.add_node(np.array([x + 0.0001, 0.0001, 0.01]))

    for i in range(node_count-1):
        beam_system.add_beam(arm_beam, i, i+1, np.array([0, 0, 1]))

    beam_system.add_force(np.array([0,0,0.5]), node_count-1)
    beam_system.add_boundary_condition(0, "en castre", 0)

    def trans_force(t):
        return np.array([0, 0, (t * (t < 0.5) + (1 - t) * (t > 0.5)) * (t < 1), 0, 0, 0])

    beam_system.solve_FEA()
    beam_system.solve_static()
    beam_system.initialize_dynamic_transient(0.02, timesteps, 0.25,0.5)
    beam_system.add_dynamic_transient_force(trans_force, node_count-1)
    beam_system.solve_dynamic_transient()

    print("F(t) = ")
    print(beam_system.nodes[node_count-1].dynamic_transient_force)

    print("U_tip(t) = ")
    print(beam_system.d_ts[:, 2 + 6*(node_count - 1)])

    print("U_d_tip(t) = ")
    print(beam_system.d_d_ts[:, 2 + 6*(node_count - 1)])

    print("Stress(t) = ")
    print(beam_system.beams[0].bending_stresses_0_ts)

    print("U_static = ")
    print(beam_system.z_displacements[node_count-1])

    print("U_tip(0.5s) = ")
    print(beam_system.d_ts[int(0.5/0.02), 2 + 6*(node_count - 1)])



if __name__ == '__main__':
    beam_30_deg()
    # timoshenko_simply_supp_point()
    # timoshenko_simply_supp()
    # timoshenko_simply_supp()
    # timoshenko_simply_supp_point_shear()
    # timoshenko_stress()
    # timoshenko_stress_circle()
    # timoshenko_freq()
    transient_displacement()
