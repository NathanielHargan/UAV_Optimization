import matplotlib.pyplot as plt
import numpy as np


def local_stiffness_3d(a, e, l, g, i_y, i_z, k, k_y=0, k_z=0):
    x = a*e/l
    phi_y = 12 * e * i_z * k_y / (a * g * l ** 2)
    y_1 = 12 * e * i_z / ((1 + phi_y) * l ** 3)
    y_2 = 6 * e * i_z / ((1 + phi_y) * l ** 2)
    y_3 = (4 + phi_y) * e * i_z / ((1 + phi_y) * l)
    y_4 = (2 - phi_y) * e * i_z / ((1 + phi_y) * l)

    phi_z = 12 * e * i_y * k_z / (a * g * l ** 2)
    z_1 = 12 * e * i_y / ((1 + phi_z) * l ** 3)
    z_2 = 6 * e * i_y / ((1 + phi_z) * l ** 2)
    z_3 = (4 + phi_z) * e * i_y / ((1 + phi_z) * l)
    z_4 = (2 - phi_z) * e * i_y / ((1 + phi_z) * l)

    s = g*k/l
    k11 = np.array([[x, 0, 0, 0, 0, 0],
                    [0, y_1, 0, 0, 0, y_2],
                    [0, 0, z_1, 0, -z_2, 0],
                    [0, 0, 0, s, 0, 0],
                    [0, 0, -z_2, 0, z_3, 0],
                    [0, y_2, 0, 0, 0, y_3]])

    k12 = np.array([[-x, 0, 0, 0, 0, 0],
                    [0, -y_1, 0, 0, 0, y_2],
                    [0, 0, -z_1, 0, -z_2, 0],
                    [0, 0, 0, -s, 0, 0],
                    [0, 0, z_2, 0, z_4, 0],
                    [0, -y_2, 0, 0, 0, y_4]])

    k21 = np.transpose(k12)

    k22 = np.array([[x, 0, 0, 0, 0, 0],
                    [0, y_1, 0, 0, 0, -y_2],
                    [0, 0, z_1, 0, z_2, 0],
                    [0, 0, 0, s, 0, 0],
                    [0, 0, z_2, 0, z_3, 0],
                    [0, -y_2, 0, 0, 0, y_3]])

    # combines the arrays
    k1 = np.concatenate([k11, k12], axis=0)
    k2 = np.concatenate([k21, k22], axis=0)
    k = np.concatenate([k1, k2], axis=1)

    return k


def transformation_matrix(coord_dir, k_node_dir):
    ortho_dir = np.cross(coord_dir, k_node_dir)

    coord_unit = coord_dir / np.linalg.norm(coord_dir)
    k_node_unit = k_node_dir / np.linalg.norm(k_node_dir)
    ortho_unit = ortho_dir / np.linalg.norm(ortho_dir)

    l1 = coord_unit[0]
    l2 = k_node_unit[0]
    l3 = ortho_unit[0]
    m1 = coord_unit[1]
    m2 = k_node_unit[1]
    m3 = ortho_unit[1]
    n1 = coord_unit[2]
    n2 = k_node_unit[2]
    n3 = ortho_unit[2]

    t = np.array([[l1, m1, n1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [l2, m2, n2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [l3, m3, n3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, l1, m1, n1, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, l2, m2, n2, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, l3, m3, n3, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, l1, m1, n1, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, l2, m2, n2, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, l3, m3, n3, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, l1, m1, n1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, l2, m2, n2],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, l3, m3, n3]])

    return t


def assemble_stiffness_3d(beam_node_indexes, ks):
    # beam_node is a nx2 matrix
    # pairs of node indexes for each element
    # [(start node,end node)
    # (start node,end node)
    # ...
    # (start node,end node)]

    # ks are a 3d array of stiffness matrixes
    # one matrix for each element i.e
    # [matrix_0,matrix_1,...,matrix_n]

    element_num = np.shape(ks)[0]

    # break it into segments
    k11s = ks[:, 0:6, 0:6]
    k12s = ks[:, 0:6, 5:-1]
    k21s = ks[:, 5:-1, 0:6]
    k22s = ks[:, 5:-1, 5:-1]

    k = np.zeros((element_num*6, element_num*6))


    for i in range(element_num): # cycles through each element

        # the stiffness matrix partitions for each element
        k11 = k11s[i]
        k12 = k12s[i]
        k21 = k21s[i]
        k22 = k22s[i]

        # the index of node 0 in the global matrix
        start_index = int(beam_node_indexes[i][0] * 6)
        # the index of node 1 in the global matrix
        end_index = int(beam_node_indexes[i][1] * 6)

        # insert element matrices into appropriate locations
        k[start_index:start_index+6, start_index:start_index+6] += k11
        k[start_index:start_index+6, end_index:end_index+6] += k12
        k[end_index:end_index+6, start_index:start_index+6] += k21
        k[end_index:end_index+6, end_index:end_index+6] += k22

    return k  # returns the global stiffness matrix


# Splits k into kuu kup kpu and kpp
def partition_stiffness_matrix(k, boundary_indexes):
    # k is the global stiffness matrix

    # boundary conditions index list
    # which index on k do boundary conditions apply

    # append new rows on k
    for i in boundary_indexes:
        k = np.concatenate((k, k[i]), axis=0)

    # delete old rows
    k = np.delete(k, boundary_indexes, axis=0)
    print("k: ", np.shape(k))
    boundary_num = len(boundary_indexes)  # number of bcs
    free_num = np.shape(k)[0] - boundary_num  # number of non-bcs
    kuu = k[0:free_num, 0:free_num]
    kup = k[0:free_num, free_num-1:-1]
    kpu = k[free_num-1:-1, 0:free_num]
    kpp = k[free_num-1:-1, free_num-1:-1]
    print("kuu: ", np.shape(kuu))
    print("kpp: ", np.shape(kpp))

    return kuu, kup, kpu, kpp


def local_to_global_stiffness_matrix(k,t):
    return np.transpose(t) @ k @ t


def solve_3D(ru, dp, kuu, kup, kpu, kpp):
    kuu_inv = np.linalg.inv(kuu)

    du = kuu_inv @ (ru - kup @ dp)
    rp = kpu @ du + kpp @ dp

    return du, rp


def shape_function(x, length):
    xi = x / length
    return np.array([1 - (3 * xi ** 2) + (2 * xi ** 3),  # N1
                     x * (1 - (2 * xi) + (xi ** 2)),  # N2
                     (3 * (xi ** 2)) - (2 * (xi ** 3)),  # N3
                     x * ((xi ** 2) - xi)])  # N4


def shape_function_derivative_1(x, length):
    xi = x / length
    return np.array([(1 / length) * (-6 * xi + 6 * xi ** 2),  # N1
                     (1 - 4 * xi + 3 * xi ** 2),  # N2
                     (1 / length) * (+6 * xi - 6 * xi ** 2),  # N3
                     x * (xi ** 2 - xi)])  # N4

