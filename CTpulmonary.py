#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import os
os.chdir('/home/szafr/ClearMap2/')
from ClearMap.Environment import * 
get_ipython().run_line_magic('gui', 'qt')
from datetime import date
import shutil


# In[ ]:


# INITIALIZE WORKSPACE
#directories and files
path = '/home/szafr/TubeMap_test/'
directory = str(path)
#filename for input
#need to run this cell twice, once with the tif file and convert and once with the npy file
filename = 'AV21_5_10um_inp.npy'
expression_raw  = filename

expression_arteries = filename

expression_auto     = filename

resources_directory = settings.resources_path

ws1 = wsp.Workspace('TubeMap', directory=directory);
ws1.update(raw=expression_raw, arteries=expression_arteries, autofluorescence=expression_auto)
ws1.info()

ws1.filename('raw')

ws1.file_list('raw')

s = ws1.source('raw')

# file conversion if we are starting with a tif

io.convert_files(ws1.file_list('raw', extension='tif'), extension='npy',
                 processes=12, verbose=True);


# In[ ]:


# ARTERIAL BINARIZATION 1

source = ws1.filename('raw');
#threshold = 450 # 750 is recommended value for vasculature
sink   = ws1.filename('binary', postfix='arteries1');
io.delete_file(sink)

binarization_parameter = vasc.default_binarization_parameter.copy();
binarization_parameter['clip']['clip_range'] = (8000,60000) #set based on min and max pixel values in the image (5000, 60000)
# (400,60000) is recommended for vasculature

binarization_parameter['lightsheet'] = None; # 0.25 is recommended value for vasculature
#binarization_parameter['deconvolve']['threshold'] = 10000 #2000
binarization_parameter['deconvolve'] = None #changing for testing image
binarization_parameter['equalize'] = None;
binarization_parameter['vesselize'] = None;
binarization_parameter['median']['selem'] = ((9,)*3) #(9, 9, 9)
#binarization_parameter['adaptive']['selem'] = (150,150,3)

processing_parameter = vasc.default_binarization_processing_parameter.copy();
processing_parameter.update(processes = 4,
                            as_memory = True, verbose=True);

vasc.binarize(source, sink,
              binarization_parameter=binarization_parameter,
              processing_parameter=processing_parameter);


# save binary image here as tif file 
io.convert_files(sink, extension='tif', processes=4, verbose=True);


# In[ ]:


#ARTERIAL BINARIZATION 1 continued

source = ws1.filename('binary', postfix='arteries1');
sink   = ws1.filename('binary', postfix='arteries1_postprocessed');
#sink_smooth = ws1.filename('binary', postfix='arteries_smoothed');

postprocessing_parameter = vasc.default_postprocessing_parameter.copy();
#postprocessing_parameter['smooth'] = dict(iterations=6) #default is 6
#postprocessing_parameter['fill'] = True #default is 6

postprocessing_processing_parameter = vasc.default_postprocessing_processing_parameter.copy();
postprocessing_processing_parameter.update(size_max = 50);

vasc.postprocess(source, sink, postprocessing_parameter=postprocessing_parameter,
                 processing_parameter=postprocessing_processing_parameter,
                 processes=None, verbose=True)

# save postprocessed binary image here
io.convert_files(sink, extension='tif', processes=None, verbose=True);

# In[ ]:


# ARTERIAL BINARIZATION 2

source = ws1.filename('raw');
#threshold = 450 # 750 is recommended value for vasculature
sink   = ws1.filename('binary', postfix='arteries2');
io.delete_file(sink)

binarization_parameter = vasc.default_binarization_parameter.copy();
binarization_parameter['clip']['clip_range'] = (3000,4000) #set based on min and max pixel values in the image
# (400,60000) is recommended for vasculature

binarization_parameter['lightsheet'] = None; # 0.25 is recommended value for vasculature
#binarization_parameter['deconvolve']['threshold'] = 2000
binarization_parameter['deconvolve'] = None #changing for testing image
binarization_parameter['equalize'] = None;
binarization_parameter['vesselize'] = None;
binarization_parameter['median'] = None
#binarization_parameter['median']['selem'] = (27, 27, 27)
#binarization_parameter['adaptive']['selem'] = (150,150,3)

processing_parameter = vasc.default_binarization_processing_parameter.copy();
processing_parameter.update(processes = 8,
                            as_memory = True, verbose=True);

vasc.binarize(source, sink,
              binarization_parameter=binarization_parameter,
              processing_parameter=processing_parameter);


# save binary image here as tif file 
io.convert_files(sink, extension='tif', processes=12, verbose=True);


# In[ ]:


#ARTERIAL BINARIZATION 2 continued

source = ws1.filename('binary', postfix='arteries2');
sink   = ws1.filename('binary', postfix='arteries2_postprocessed');
#sink_smooth = ws1.filename('binary', postfix='arteries_smoothed');

postprocessing_parameter = vasc.default_postprocessing_parameter.copy();
#postprocessing_parameter['smooth'] = dict(iterations=6) #default is 6

postprocessing_processing_parameter = vasc.default_postprocessing_processing_parameter.copy();
postprocessing_processing_parameter.update(size_max = 50);

vasc.postprocess(source, sink, postprocessing_parameter=postprocessing_parameter,
                 processing_parameter=postprocessing_processing_parameter,
                 processes=None, verbose=True)

# save postprocessed binary image here
io.convert_files(sink, extension='tif', processes=None, verbose=True);
#io.convert_files(sink, extension='tif', processes=None, verbose=True);

#%% Combine binaries
  
source          = ws.filename('binary', postfix='arteries1_postprocessed');
source_arteries = ws.filename('binary', postfix='arteries2_postprocessed');
sink            = ws.filename('binary', postfix='final');
  
bp.process(np.logical_or, [source, source_arteries], sink, size_max=500, overlap=0, processes=None, verbose=True)
  
#p2d.plot([source, source_arteries, sink]);
  
  

# In[ ]:
io.convert_files(path + 'binary_arteries1_postprocessed_fiji.tif', extension='npy',
                 processes=12, verbose=True);
fiji_file = np.load(path + 'binary_arteries1_postprocessed_fiji.npy')
fiji_file = fiji_file > 0
np.save(path + 'binary_arteries1_postprocessed_fiji.npy', fiji_file)
del fiji_file

#%%
binary   = ws1.filename('binary', postfix='arteries1_postprocessed');
skeleton = ws1.filename('skeleton')

skl.skeletonize(binary, sink=skeleton, delete_border=True, verbose=True);

graph_raw = gp.graph_from_skeleton(ws1.filename('skeleton'), verbose=True)
#graph_raw.save(ws.filename('graph', postfix='raw'))

#ws1.filename('raw')
#OR binary
coordinates = graph_raw.vertex_coordinates();
radii, indices = mr.measure_radius(binary, coordinates,
    #                               value=0, fraction=None, max_radius=150,
                                   value=None, fraction=0.10, max_radius=100,
                                   return_indices=True, default=-1, verbose=True);
graph_raw.set_vertex_radii(radii)

graph_raw.save(ws1.filename('graph', postfix='raw'))

# Graph cleaning
graph_cleaned = gp.clean_graph(graph_raw,
                               vertex_mappings = {'coordinates'   : gp.mean_vertex_coordinates,
                                                  'radii'         : np.max,
                                                  'artery_binary' : np.max,
                                                  'artery_raw'    : np.max},
                               verbose=True)

# Save cleaned graph

graph_cleaned.save(ws1.filename('graph', postfix='cleaned'))
#graph_cleaned = grp.load(ws.filename('graph', postfix='cleaned'));

# Graph reduction

def vote(expression):
  return np.sum(expression) >= len(expression) / 1.5;

graph_reduced = gp.reduce_graph(graph_cleaned, edge_length=True,
                          edge_to_edge_mappings = {'length' : np.sum},
                          vertex_to_edge_mappings={'artery_binary' : vote,
                                                   'artery_raw'    : np.max,
                                                   'radii'         : np.max},
                          edge_geometry_vertex_properties=['coordinates', 'radii', 'artery_binary', 'artery_raw'],
                          edge_geometry_edge_properties=None,
                          return_maps=False, verbose=True)

graph_reduced.save(ws1.filename('graph', postfix='reduced'))
#graph_reduced = grp.load(ws.filename('graph', postfix='reduced'));

# Graph largests component

graph = graph_reduced.largest_component()
graph.save(ws1.filename('graph', postfix='reduced_largest'))

graph_reduced_sk = gp.graph_to_skeleton(graph)
np.save(path + 'reduced_skeleton.npy', graph_reduced_sk)
io.convert_files(path + 'reduced_skeleton.npy', extension='tif', processes=None, verbose=True);

graph_reduced = grp.load(ws1.filename('graph', postfix='reduced_largest'));
verticies = graph_reduced.vertex_coordinates()
connectivity = graph_reduced.edge_connectivity()
radii = graph_reduced.vertex_radii()

print(connectivity)
print(len(connectivity))
print(verticies)
print(radii)
print(len(radii))

conn_save = directory +  "connectivity_reduced.p"
vert_save = directory +  "verticies_reduced.p"
rad_save = directory +  "radii_reduced.p"

import pickle
pickle.dump(connectivity, open( conn_save, "wb" ) )
pickle.dump(verticies, open( vert_save, "wb" ) )
pickle.dump(radii, open( rad_save, "wb" ) )

##%%
##import ClearMap.Analysis.Graphs.GraphProcessing as gp
#import ClearMap.Analysis.Graphs.GraphGt as ggt
##import ClearMap.Analysis.Graphs.Graph as grp
#g_test = ggt.Graph()
#needle_vertices = [11602, 16068]
##g = gp.ggt.Graph(graph_reduced)
#graph_reduced.remove_vertex(needle_vertices)
#g_test = graph_reduced.largest_component()
##graph_trim2 = gp.reduce_graph(g_trim1, edge_length=True,
##                          edge_to_edge_mappings = {'length' : np.sum},
##                          vertex_to_edge_mappings={'artery_binary' : vote,
##                                                   'artery_raw'    : np.max,
##                                                   'radii'         : np.max},
##                          edge_geometry_vertex_properties=['coordinates', 'radii', 'artery_binary', 'artery_raw'],
##                          edge_geometry_edge_properties=None,
##                          return_maps=False, verbose=True)
##
##g_trim2 = g_trim1.largest_component()
#
#verticies = g_test.vertex_coordinates()
#connectivity = g_test.edge_connectivity()
#radii = g_test.vertex_radii()
#
#pickle.dump(connectivity, open( "/home/jszafron/Documents/YaleHypoxia/TubeMap_test/connectivity_new.p", "wb" ) )
#pickle.dump(verticies, open( "/home/jszafron/Documents/YaleHypoxia/TubeMap_test/verticies_new.p", "wb" ) )
#pickle.dump(radii, open( "/home/jszafron/Documents/YaleHypoxia/TubeMap_test/radii_new.p", "wb" ) )

