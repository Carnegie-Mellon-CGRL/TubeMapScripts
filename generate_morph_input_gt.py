#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 07:48:29 2022

@author: jszafron
"""
from datetime import date
from graph_tool.all import *
import shutil
import numpy as np
import math
import random
import pickle
import reduce_graph_helper
import matplotlib.pyplot as plt

path = '/home/jszafron/Documents/source/PAMorphometry/AVShunt4wk/AV21_5_10um/TubeMap_test/'
directory = str(path)

import os
import sys
file_dir = os.path.dirname(directory)
sys.path.append(file_dir)

import graphToCenterline

# Load inital graph from segmentation code
#path = '/home/jszafron/Documents/source/PAMorphometry/AV30_4_10um/TubeMap_test/'
#directory = str(path)

ModelName = 'AV21_5'

#%%
g_init = load_graph(path + 'graph_reduced_largest.gt')

inlet_vert1 = 60
inlet_vert2 = 63
#inlet_vert3 = 13758
g_init.vp.radii[inlet_vert1] = 52
g_init.vp.radii[inlet_vert2] = 52
#g_init.vp.radii[inlet_vert3] = 52 

remove_self_loops(g_init)
remove_parallel_edges(g_init)

#clean up the input graph
#g_trim = g_init #for needle extraction

# Remove segments resulting from needle insertion to MPA
# Requires manual identification of needle nodes

#Set upstream MPA node to same diameter as downstream node
# ds_radius = g_init.vp.radii[11429]
# g_init.vp.radii[16035] = ds_radius

# #Keep only the right PA tree
# ds_radius = g_init.vp.radii[15621]
# g_init.vp.radii[11429] = ds_radius

#Keep only the left PA tree - No artifact in LPA
# ds_radius = g_init.vp.radii[15621]
# g_init.vp.radii[11429] = ds_radius

g_trim = Graph()
remove_vertices = [81, 13758, 86, 87, 88, 145, 225, 394, 413, 435, 533,
                   144, 145, 199, 233, 246, 296, 318, 357, 439, 13842,
                   226, 466, 473, 530, 615, 646, 730, 774, 820, 865, 911,
                   912, 980, 1108, 1372, 13889, 13972, 255, 1115, 13944] #for MPA

##268, 286, 327, 362, 431, 11991, 12023, 12122 
g_init.remove_vertex(remove_vertices)
g_trim = extract_largest_component(g_init, None, True)

#remove_vertices2 = [50, 66, 108, 119, 420, 450, 591, 720, 735, 8620, 8641]
#g_trim.remove_vertex(remove_vertices2)
#g_trim = extract_largest_component(g_trim, None, True)

# #Optional code to remove 
# #remove vertexes smaller than threshold
# min_vert_size = 3 # in px, ~0.03 mm / pix
# remove_ind = []
# for v in g_trim.iter_vertices():
    # if g_trim.vp.radii[v] < min_vert_size:
        # remove_ind.append(v)

# g_trim.remove_vertex(remove_ind)
# g_trim = extract_largest_component(g_trim, None, True)

#find new input vertex ind
curr_max_rad = 0
curr_max_ind = 0
for v in g_trim.iter_vertices():
    if g_trim.get_total_degrees([v]) == 1:
        if g_trim.vp.radii[v] > curr_max_rad:
            curr_max_rad = g_trim.vp.radii[v]
            curr_max_ind = v
input_vertex_ind = curr_max_ind
in_vertex = g_trim.vertex(input_vertex_ind)

#Remove small branches near the inlet
#Loop through all vertices, find order 1 vertices
#Find distance to inlet < l_path_thres
#store vertex ID to be removed
#remove vertices
#reduce the graph
max_dist = 15
vert_ind_prune = []
for v in g_trim.iter_vertices():
    if g_trim.get_total_degrees([v]) == 1:
        path2inlet = shortest_path(g_trim, v, in_vertex)
        if len(path2inlet[0]) < max_dist and v != in_vertex:
            vert_ind_prune.append(v)
            
g_trim.remove_vertex(vert_ind_prune)

#iterate to check if there are more large branches to clean
count = 0
max_count = 4
bool_large_branches = 1
while bool_large_branches and count < max_count:
    #max_dist = math.floor(max_dist / 2)
    #find new input vertex ind
    curr_max_rad = 0
    curr_max_ind = 0
    for v in g_trim.iter_vertices():
        if g_trim.get_total_degrees([v]) == 1:
            if g_trim.vp.radii[v] > curr_max_rad:
                curr_max_rad = g_trim.vp.radii[v]
                curr_max_ind = v
    input_vertex_ind = curr_max_ind
    in_vertex = g_trim.vertex(input_vertex_ind)

    #Remove small branches near the inlet
    vert_ind_prune = []
    for v in g_trim.iter_vertices():
        if g_trim.get_total_degrees([v]) == 1:
            path2inlet = shortest_path(g_trim, v, in_vertex)
            if len(path2inlet[0]) < max_dist and v != in_vertex:
                vert_ind_prune.append(v)
    if len(vert_ind_prune) > 0:            
        g_trim.remove_vertex(vert_ind_prune)
    else:
        bool_large_branches = 0
        

##Optional code to remove 
##remove edges from vertexes with degree > 3
#max_seg_per_node = 3
#for v in g_trim.iter_vertices():
#    if len(g_trim.get_all_edges(v)) > max_seg_per_node:
#        curr_vert = g_trim.vertex(v)
#        iter_n = len(g_trim.get_all_edges(v)) - max_seg_per_node
#        #remove the edges to the smallest vertices
#        for i in range(iter_n):
#            smallest_rad = 10000
#            smallest_rad_ind = 0
#            for v_n in curr_vert.all_neighbors():
#                neighbor_ind = g_trim.vertex_index[v_n]
#                if g_trim.vp.radii[neighbor_ind] < smallest_rad:
#                    smallest_rad = g_trim.vp.radii[neighbor_ind]
#                    smallest_rad_ind = neighbor_ind
#            remove_edge = g_trim.edge(v, smallest_rad_ind)
#            g_trim.remove_edge(remove_edge)

g_trim.reindex_edges()
g_trim = extract_largest_component(g_trim, None, True)

edge_rad = g_trim.ep.radii
edge_rad_arr = edge_rad.get_array()
edge_rad_arr_inv = 1 / np.array(edge_rad_arr)
edge_prop_rad_arr_inv = g_trim.new_edge_property('double', vals = edge_rad_arr_inv)
g_trim.ep["inv_rad"] = edge_prop_rad_arr_inv
min_tree_edges = min_spanning_tree(g_trim, weights = g_trim.ep["inv_rad"], root = in_vertex)

g_trim.set_edge_filter(min_tree_edges)
g_trim.purge_edges()
g_trim = extract_largest_component(g_trim, None, True)

vertex_coords = []
radii = []
conn = []
for v in g_trim.iter_vertices():
    vertex_coords.append(g_trim.vp.coordinates[v])
    radii.append(g_trim.vp.radii[v])

conn_save = directory +  "connectivity_pruned.p"
vert_save = directory +  "verticies_pruned.p"
rad_save = directory +  "radii_pruned.p"

conn = g_trim.get_edges()

pickle.dump(conn, open( conn_save, "wb" ) )
pickle.dump(vertex_coords, open( vert_save, "wb" ) )
pickle.dump(radii, open( rad_save, "wb" ) )

g_trim.save("graph_pruned.gt", fmt="gt")

#Reduce graph
reduce_graph_helper.main(path)

#Generate intermediate vtk graph
graphToCenterline.main(path, '_reduced_pruned')

#Load reduced graph
g_prune = load_graph(path + 'graph_pruned_reduced.gt')

#%%
#find new input vertex ind
curr_max_rad = 0
curr_max_ind = 0
for v in g_prune.iter_vertices():
    if g_prune.get_total_degrees([v]) == 1:
        if g_prune.vp.radii[v] > curr_max_rad:
            curr_max_rad = g_prune.vp.radii[v]
            curr_max_ind = v
input_vertex_ind = curr_max_ind
in_vertex = g_prune.vertex(input_vertex_ind)

#Find outlets
out_vert_inds = []
for v in g_prune.iter_vertices():
    if g_prune.get_total_degrees([v]) == 1:
        out_vert_inds.append(v)            

n_paths = 200
random_out = random.sample(out_vert_inds, n_paths)

#curr_vert = g_prune.vertex(random_out[0])
#path2inlet = shortest_path(g_prune, curr_vert, in_vertex)

#%%
n_bins = 5

path2inlet = []
radius_gen = [[]]
radius_bins = [[0] * n_bins]
fig, ax = plt.subplots(1,2)
for i in range(n_paths):
    print('Branch: ', i)
    curr_vert = g_prune.vertex(random_out[i])
    curr_path = shortest_path(g_prune, curr_vert, in_vertex)
    
    #if len(curr_path > 2 * )
    path2inlet.append(curr_path)
    
    x = range(len(path2inlet[i][0][:]))
    for j in x:
        curr_rad = g_prune.vp.radii[path2inlet[i][0][j]]
        radius_gen[i].append(curr_rad)
    
    y = radius_gen[i][:]
    ax[0].plot(x,y)
    
    #create normalized branch with n_bins
    branch_bin_size = math.floor(len(radius_gen[i]) / n_bins)
    for k in range(n_bins - 1):
        print('Bin: ', k, ' Range: ', k * branch_bin_size, ' to ', (k + 1) * branch_bin_size - 1 * (branch_bin_size > 1))
        radius_bins[i][k] = np.mean(radius_gen[i][k * branch_bin_size: (k + 1) * branch_bin_size - 1 * (branch_bin_size > 1)])
    #remaining vertices into last bin
    print('Bin: ', n_bins - 1, ' Range: ', branch_bin_size * (n_bins - 1), ' to ', len(radius_gen[i][:]))
    radius_bins[i][n_bins - 1] = np.mean(radius_gen[i][branch_bin_size * (n_bins - 1):])
    
    if i < n_paths - 1:     
        radius_gen.append([])
        radius_bins.append([0] * n_bins)
    print('-------------')

arr_radius_gen = np.array(radius_gen)
arr_radius_bins = np.array(radius_bins)
mean_radius_bins = np.average(arr_radius_bins, axis = 0)
std_radius_bins = np.std(arr_radius_bins, axis = 0)

ax[1].plot(range(n_bins), mean_radius_bins)
plt.show()

np.save(path + ModelName + '_radius_gens.npy', arr_radius_gen)
np.save(path + ModelName + '_radius_bins.npy', arr_radius_bins)

#%%
#Create Strahler input file
# Open file
file = open(path + ModelName + ".in", "w")
file2 = open(path + "groupID_avgArea.dat", "w")

# Write header
file.write("# ================================\n# " + ModelName + " MODEL - UNITS IN CGS\n# ================================\n\n")

#img resolution 10 um/pix
voxel_2_mm = 10 #mm/vx scaling for uCT image
mm_2_cm = 1 #conversion cm/mm

# Node Header
file.write("\n\n### DO NOT CHANGE THIS SECTION - generated automatically \n#")
file.write("\n# ==========\n# NODE CARD\n# ==========\n# - Node Name (double)\n# - Node X Coordinate (double)\n# - Node Y Coordinate (double)\n# - Node Z Coordinate (double)\n\n")

#Write out the node values
for v in g_prune.vertices():
    xyz = g_prune.vp.coordinates[v]
    file.write("NODE " + str(v) + " " + str(xyz[0] * voxel_2_mm * mm_2_cm) + " " + str(xyz[1]* voxel_2_mm * mm_2_cm) + " " + str(xyz[2] * voxel_2_mm * mm_2_cm) + "\n")

# Joint Header
file.write("\n\n\n### DO NOT CHANGE THIS SECTION - generated automatically \n#")
file.write("\n# ==========\n# JOINT CARD\n# ==========\n# - Joint Name (string)\n# - Joint Node (double)\n# - Joint Inlet Name (string)\n# - Joint Outlet Name (string)\n\n")
# JointInlet and JointOutlet Header
file.write("\n### DO NOT CHANGE THIS SECTION - generated automatically \n#")
file.write("\n# ================================\n# JOINTINLET AND JOINTOUTLET CARDS\n# ================================\n# - Inlet/Outlet Name (string)\n# - Total Number of segments (int)\n# - List of segments (list of int)\n\n")

#Write out the joint values
joint_count = 0;
distances = shortest_distance(g_prune, in_vertex)
distance_arr = distances.get_array()
for v in g_prune.iter_vertices():
    if g_prune.get_total_degrees([v]) > 1:
        
        joint_count += 1 #increment the joint #
        nb_vert = g_prune.get_all_neighbors(v)
        
        #Set the inlet to the segment with the shortest path
        inlet_vert = nb_vert[np.argmin(distance_arr[nb_vert])]
        inlet_edge = g_prune.edge(inlet_vert, v)
        inlet_edge_ind = g_prune.edge_index[inlet_edge]
        
        #Gather the outlet segments
        outlet_vert = np.delete(nb_vert, np.argmin(distance_arr[nb_vert]))
        num_outlet = len(outlet_vert)
        outlet_edge = []
        outlet_edge_ind = []
        for ov in range(num_outlet):
            curr_outlet_edge = g_prune.edge(v, outlet_vert[ov])
            outlet_edge.append(curr_outlet_edge)
            curr_outlet_edge_ind = g_prune.edge_index[curr_outlet_edge]
            outlet_edge_ind.append(curr_outlet_edge_ind)
        
        file.write("JOINT J" + str(joint_count) + " " + str(v) + " IN" + str(joint_count) + " OUT" + str(joint_count) + "\n")
    
        file.write("JOINTINLET IN" + str(joint_count) + " " + "1 " + str(inlet_edge_ind) + "\n")
        file.write("JOINTOUTLET OUT" + str(joint_count) + " " + str(num_outlet))
        
        for j in range(num_outlet): 
            file.write(" " + str(outlet_edge_ind[j]))
        
        file.write("\n\n")
        
# Segment Header
file.write("# ============\n# SEGMENT CARD\n# ============\n# - Segment Name (string)\n# - Segment ID (int)\n# - Segment Length (double)\n# - Total Finite Elements in Segment (int)\n# - Segment Inlet Node (int)\n# - Segment Outlet Node (int)\n# - Segment Inlet Area (double)\n# - Segment Outlet Area (double)\n# - Segment Inflow Value (double)\n# - Segment Material (string)\n# - Type of Loss (string - 'NONE')\n# - 0.0 (double)\n# - 0 (int)\n# - 0 (int)\n# - Boundary Condition Type (string - 'NOBOUND','RESISTANCE','RCR')\n# - Data Table Name (string)\n\n")  
edge_count = 0
for e in g_prune.iter_edges():

    matname = "MAT1"
    numfe = 0
    
    curr_edge = g_prune.edge(e[0], e[1])
    edge_ind = g_prune.edge_index[curr_edge]
    
    #Check if inlet and print info for quantify morphometry
    if (e[0] == in_vertex or e[1] == in_vertex):
        print("Group"+ str(edge_count))
        print(g_prune.vp.radii[e[0]])
        print(g_prune.vp.radii[e[1]])

    coord_v1 = g_prune.vp.coordinates[e[0]]
    coord_v2 = g_prune.vp.coordinates[e[1]]
    length = ( (coord_v1[0] - coord_v2[0])**2 + (coord_v1[1] - coord_v2[1])**2 + (coord_v1[2] - coord_v2[2])**2 )**(1 / 2) * voxel_2_mm * mm_2_cm
    
    #mean of nodal radii
    rn1 = g_prune.vp.radii[e[0]]
    Ai = math.pi * (rn1 * voxel_2_mm * mm_2_cm)**2
    rn2 = g_prune.vp.radii[e[1]]
    Ao = math.pi * (rn2 * voxel_2_mm * mm_2_cm)**2
    Area_avg = (Ai + Ao) / 2
    Area_std = np.std([Ai, Ao])
    
    file.write("SEGMENT" + " " + "Group"+ str(edge_count)+"_Seg"+str(edge_ind) + " " + str(edge_ind) + " "+ str(length) + " " + str(numfe) + " "+ str(e[0]) + " " + str(e[1]) + " " + str(Ai)+ " " + str(Ao)+ " " +"0.0 "+ matname + " NONE 0.0 0 0 NOBOUND NONE \n")
    file2.write("Group"+ str(edge_count)+"_Seg"+str(edge_ind) + " " + str(Area_avg) + " " + str(Area_std) + "\n")
    edge_count += 1

file.close()
file2.close()