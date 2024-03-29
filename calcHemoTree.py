#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 12:22:36 2024

@author: szafr
"""

from datetime import date
from graph_tool.all import *
import shutil
import numpy as np
import math
import random
import pickle
#import pandas as pd
#import reduce_graph_helper
import matplotlib.pyplot as plt

def main(Path, ModelName, Q, P):
    g_prune = load_graph(Path + ModelName + '_graph_pruned_reduced.gt')
    #img resolution 10 um/pix
    voxel_2_mm = 10 #mm/vx scaling for uCT image
    mm_2_cm = 1 #conversion cm/mm
    
    #vars
    mu = 0.04 #dynamic viscosity of blood
    
    
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
    
    edge_seg_res = []
    edge_count = 0       
    for e in g_prune.iter_edges():
        
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
        rn2 = g_prune.vp.radii[e[1]]
        r_mean = (rn2 + rn1)/2
        
        #calc seg res
        seg_res = 8 * mu * length / (np.pi * r_mean**4)
        edge_seg_res.append(seg_res)
        
        edge_count += 1
    
    
if __name__ == '__main__':
    # executed as script
    Q = 1
    P = 1
    Path = '/home/szafr/AVResults/AVResults/ResultsCompare/'
    ModelName = 'Shunt4wk_AV21_5'
    main(Path, ModelName, Q, P)