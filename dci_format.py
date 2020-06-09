__author__ = 'Mathis Messager'

import arcpy
import os
import glob
import itertools
from collections import defaultdict

#Set analysis parameters
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = True
arcpy.env.XYResolution = "0.001 meters"
arcpy.env.XYTolerance = "0.001 meters"

#Input data and variables
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
blevels = ['04', '06', '08', '10']
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')
crsSIRGAS = arcpy.SpatialReference(5880) #SIRGAS_2000_Brazil_Polyconic (ESPG: 5880)
streamatlas = os.path.join(datadir, 'HydroATLAS/RiverATLAS_v10.gdb/RiverATLAS_v10')
br_admin = os.path.join(datadir, 'ibge/BR/BRUFE250GC_SIR.shp')

#Output variables
ngacoast = os.path.join(resdir, "ngacoast.shp")
br_bound = os.path.join(resdir, "br_boundaries.shp")
dcigdb = os.path.join(resdir, 'dci.gdb')
if not os.path.exists(dcigdb):
    arcpy.CreateFileGDB_management(resdir, 'dci')
edit = arcpy.da.Editor(dcigdb)

snbr_original = os.path.join(dcigdb, 'streamATLAS_Brazil')
netftdat = os.path.join(dcigdb, 'streamATLAS_geonetwork')
snbr_proj =  os.path.join(netftdat,'streamATLAS_Brazil_proj')
snbr_startvert = os.path.join(dcigdb, 'streamATLAS_startvert')
snbr_startsplit = os.path.join(dcigdb, 'streamATLAS_startsplit')
snbr_startsplitvert = os.path.join(dcigdb, 'streamATLAS_startsplitvert')
geonetini = os.path.join(netftdat, 'streamATLAS_geonet_ini')
snbr_trim = os.path.join(dcigdb, 'streamATLAS_Brazil_startsplit_trim')
snbr_clean = os.path.join(dcigdb,'streamATLAS_Brazil_clean')

basinbr4 = os.path.join(dcigdb, "BasinATLAS_lev04_v10_br")
netproj4 = os.path.join(netftdat, 'streamATLAS_proj_level04')
netproj10 = os.path.join(netftdat, 'streamATLAS_proj_level10')
geonet4 = os.path.join(netftdat, 'streamATLAS_geonet_level4')
geonet10 = os.path.join(netftdat, 'streamATLAS_geonet_level10')

dams_original = os.path.join(datadir, 'DamsBrazil/Aproveitamento_Hidrelétricos_AHE.shp')
dams_originalgdb = os.path.join(dcigdb, 'Aproveitamento_Hidreletricos_AHE')
damsel = os.path.join(dcigdb, 'damfull_simplified_300m')
damsbuf300 = os.path.join(dcigdb, 'damsbuf300')
damsbufinter = os.path.join(dcigdb, 'damsbufinter')
damsnap = os.path.join(dcigdb, 'damfull_snapped')
damsnap_join = os.path.join(dcigdb, 'damfull_snappedjoin')
damsnap_edit = os.path.join(dcigdb, 'damfull_snappededit')
damsnap_editclean = os.path.join(dcigdb, 'damfull_snappededitclean')
damsnap_editcleanjoin = os.path.join(dcigdb, 'damfull_join')

outsnsplit = os.path.join(netftdat, 'streamATLAS_split')
geonetout = os.path.join(netftdat,'streamATLAS_geonet')

netp = os.path.join(dcigdb, 'netpointstart')
netendp = os.path.join(dcigdb, 'netpointend')
netdanglep = os.path.join(dcigdb, 'netpointdangle')
net_diss = os.path.join(dcigdb, 'netproj_diss')

dam_tab = os.path.join(dcigdb, 'damattributes')
netseg_tab = os.path.join(dcigdb, 'networkattributes')

########################################################################################################################
# GENERAL DATA PRE-FORMATTING
########################################################################################################################
#Pre-format coastline vector
arcpy.Merge_management(glob.glob(os.path.join(datadir, 'NGA/*/cd*.shp')), ngacoast)

#Dissolve Brazil municipalities
arcpy.Dissolve_management(br_admin, br_bound)

for level in blevels:
    #---- Filter basins and rivers of Brazil ----
    basin = os.path.join(datadir, "HydroATLAS/BasinATLAS_v10.gdb/BasinATLAS_lev{}_v10".format(level))
    basinbr = os.path.join(dcigdb, "{}_br".format(os.path.split(basin)[1]))
    basinbrinter = os.path.join(dcigdb, "{}_brinter".format(os.path.split(basin)[1]))
    basinbrproj = basinbr + 'proj'
    if not arcpy.Exists(basinbr):
        print('Filter BasinATLAS for Brazil for level {}...'.format(level))
        if not 'AREA_GEOFULL' in [f.name for f in arcpy.ListFields(basin)]:
            arcpy.AddGeometryAttributes_management(basin, 'AREA_GEODESIC', Area_Unit = "SQUARE_KILOMETERS")
            arcpy.AlterField_management(basin, 'AREA_GEO', new_field_name='AREA_GEOFULL', new_field_alias = 'AREA_GEOFULL')
        arcpy.Intersect_analysis([basin, br_bound], basinbrinter)
        arcpy.AddGeometryAttributes_management(basinbrinter, 'AREA_GEODESIC', Area_Unit = 'SQUARE_KILOMETERS')
        #Only keep basins that have at least 1% of their area within Brazil
        arcpy.MakeFeatureLayer_management(basinbrinter, 'basinbrinterlyr',
                                          where_clause= '(AREA_GEO/AREA_GEOFULL) > 0.01')
        basinbrlist = [row[0] for row in arcpy.da.SearchCursor('basinbrinterlyr', ['HYBAS_ID'])]
        arcpy.MakeFeatureLayer_management(basin, 'basinbrlyr', where_clause= 'HYBAS_ID IN {}'.format(tuple(basinbrlist)))
        arcpy.CopyFeatures_management('basinbrlyr', basinbr)

        #Subset StreamATLAS based on level 4 basins that intersect Brazil ----
        if level == '04':
            arcpy.Clip_analysis(streamatlas, 'basinbrlyr', snbr_original)

        #Project basins
        arcpy.Project_management(basinbr, basinbrproj, crsSIRGAS)

#---- Remove suspended and non-licensed dams ----
arcpy.MakeFeatureLayer_management(dams_original, 'dams_sublyr')
arcpy.SelectLayerByAttribute_management('dams_sublyr', 'NEW_SELECTION', where_clause= "NOT {0} IN ('Desativado' , 'Revogado', 'Extinta')".format('"ESTAGIO_1"'))

#-----Clean and project river network + project dams ----
#Identify and remove the bugs in the streamATLAS network
#Project dam dataset and import it into gdb
if not arcpy.Exists(dams_originalgdb):
    arcpy.Project_management('dams_sublyr', dams_originalgdb, crsSIRGAS)

#Create feature dataset
#pr.XYTolerance = 0.01
if not arcpy.Exists(netftdat):
    arcpy.CreateFeatureDataset_management(dcigdb, os.path.split(netftdat)[1], crsSIRGAS) #pr
#arcpy.Describe(netftdat).SpatialReference.XYTolerance

#Project network and import into feature dataset
arcpy.Project_management(snbr_original, snbr_proj, crsSIRGAS)

######Remove bifurcations in network (start of topological loops)
arcpy.FeatureToLine_management(snbr_proj, snbr_startsplit)
arcpy.FeatureVerticesToPoints_management(snbr_startsplit, snbr_startsplitvert, point_location='START')

vertIDs = defaultdict(list)
with arcpy.da.SearchCursor(snbr_startsplitvert, ['Shape@XY', 'REACH_ID']) as cursor:
    x=0
    for row in cursor:
        if x%100000 == 0:
            print('Processed {} records for duplicate shapes...'.format(x))
        vertIDs[row[0]].append(row[1])
        x+=1
dupliverts = [v for k,v in vertIDs.iteritems() if len(v)>1]

#Inspect them
#When one reach's start point is in the middle of the other reach, it will split the latter in two. So shorten the one
#that has only one corresponding segment. When the two reaches share the same start point, simply shorten the first one
dupliverts_flat = [v2 for v1 in dupliverts for v2 in v1]
#Count number of splitted segments associated with each reach
dupliverts_lendict = defaultdict(int)
with arcpy.da.SearchCursor(snbr_startsplit, ['REACH_ID']) as cursor:
    for row in cursor:
        if row[0] in dupliverts_flat:
            dupliverts_lendict[row[0]] += 1

#When one of a pair of bifurcating reaches has 2 segments, select the other
trimreachlist = []
for startp in dupliverts:
    check2 = [dupliverts_lendict[p] == 2 for p in startp]
    if sum(check2) == 1:
        trimreachlist.append([i for (i, v) in zip(startp, check2) if not v][0])
    else:
        trimreachlist.append(startp[0])

#Create buffers around start vertices of selected segments
sel_expr = 'REACH_ID IN {}'.format(tuple(trimreachlist))
arcpy.MakeFeatureLayer_management(snbr_startsplitvert, 'startsplitvertlyr', where_clause=sel_expr)
arcpy.Buffer_analysis('startsplitvertlyr', os.path.join(dcigdb, 'bifurcbuf'), '50 meters', method='GEODESIC')

#Trim selected segments with buffers
arcpy.MakeFeatureLayer_management(snbr_startsplit, 'startsplitlyr')
arcpy.SelectLayerByAttribute_management('startsplitlyr', 'NEW_SELECTION', sel_expr)
arcpy.Erase_analysis('startsplitlyr', os.path.join(dcigdb, 'bifurcbuf'), os.path.join(dcigdb, 'snsubtrim'))

#Merge trimmed segments with rest of network and re-dissolve it by REACH_ID
arcpy.SelectLayerByAttribute_management('startsplitlyr', 'NEW_SELECTION',
                                        where_clause='REACH_ID NOT IN {}'.format(tuple(trimreachlist)))
arcpy.Merge_management(['startsplitlyr', os.path.join(dcigdb, 'snsubtrim')], snbr_trim)
arcpy.Dissolve_management(snbr_trim, snbr_clean, dissolve_field = 'REACH_ID')

########################################################################################################################
# BARRIER PREPARATION
########################################################################################################################
#14) Clean the dam dataset
# (“Aproveitamento_Hidrelétricos_AHE.shp” – small and large/current and future).
# Merge centrals (e.g. Santa Lucia I and Santa Lucia II) and remove old ones (e.g. Dallasta, São Domingos and Ludesa).

#Create unique dam ID
print('Create unique dam ID...')
arcpy.AddField_management(dams_originalgdb, 'DAMID', 'TEXT')
with arcpy.da.UpdateCursor(dams_originalgdb, [arcpy.Describe(dams_originalgdb).OIDFieldName, 'DAMID']) as cursor:
    for row in cursor:
        row[1] = 'D{}'.format(row[0])
        cursor.updateRow(row)

#create a buffer layer (300m) of the dams
print('Remove overlapping dams...')
buffdist = 300
arcpy.Buffer_analysis(dams_originalgdb, damsbuf300, '{} meters'.format(buffdist), line_side='FULL',
                      dissolve_option='NONE', method='GEODESIC')

# Second, identify what dam are within 300 m of each other and only keep the larger ones
arcpy.Intersect_analysis([dams_originalgdb, damsbuf300], damsbufinter, join_attributes='ALL')
arcpy.Sort_management(damsbufinter, damsbufinter + 'sort', sort_field=[['POT_KW', 'DESCENDING']])

# [f.name for f in arcpy.ListFields(damsbufinter + 'sort')]
keepdict = defaultdict(list)
outdict = []
with arcpy.da.SearchCursor(damsbufinter + 'sort', ['DAMID', 'DAMID_1']) as cursor:
    for row in cursor:
        if row[1] not in outdict: #if dam not already excluded
            if row[0] not in keepdict: #if dam buffer not already in the dictionary
                keepdict[row[0]] = row[1]
            else: #if dam buffer already has a dam associated with it (largest, given that they were sorted in descending size order)
                outdict.append(row[1])

print('Deleting {0} dams within {1} meters from a larger dam...'.format(len(outdict), buffdist))
arcpy.CopyFeatures_management(dams_originalgdb, damsel) #Seems that the FID shifted by one during copy
damkeepset = set(keepdict.values())
with arcpy.da.UpdateCursor(damsel, ['DAMID']) as cursor:
    for row in cursor:
        if row[0] not in damkeepset:
            print(row[0])
            cursor.deleteRow()

# 16) Snap the dams with the streamATLAS-debug networks
#Check distance from dams to river network
print('Snap dams to river network...')
arcpy.CopyFeatures_management(damsel, damsnap)
arcpy.Near_analysis(damsnap, snbr_clean, location = 'LOCATION', angle = 'NO_ANGLE')

#Snap
snap_env = [snbr_clean, "EDGE", "200 Meters"]
arcpy.Snap_edit(damsnap, [snap_env])

#Compare dam size with modelled discharge of stream it was snapped to
arcpy.SpatialJoin_analysis(damsnap, streamatlas, damsnap_join, join_operation= 'JOIN_ONE_TO_ONE',
                           join_type = 'KEEP_ALL', match_option = 'INTERSECT', search_radius = '10 meters')
arcpy.AddField_management(damsnap_join, 'KWTOQ', 'FLOAT')
arcpy.CalculateField_management(damsnap_join, 'KWTOQ', '!POT_KW!/!dis_m3_pyr!', expression_type= 'PYTHON')

#Compare drainage areas between the dam database and that obtained with HydroSHEDS
arcpy.AddField_management(damsnap_join, 'DRENCOMP', 'FLOAT')
arcpy.CalculateField_management(damsnap_join, 'DRENCOMP', '!AREA_DREN!/!UPLAND_SKM!', expression_type= 'PYTHON')

#Copy dams to edit manually - perform all edits on damsnap_edit
arcpy.CopyFeatures_management(damsnap_join, os.path.join(datadir, damsnap_edit))
#Check those dams 200-500 m from river network with imagery and water mask, large hydropower (Tipo_1 = UHE), and those with high/low KWTOQ.
#Either move them to the network or delete them.
#Create field to know whether dam was checked and edited (1) or not edited (0), otherwise (if not checked) 'null'
arcpy.AddField_management(damsnap_edit, 'manualsnap', 'TEXT')
arcpy.AddField_management(damsnap_edit, 'comment', 'TEXT')

#For dams in operation, the location of a dam, reservoir, or plant was the primary criterion of placement.
#In the absence of visible structure in the vicinity...etc.
#The authors estimate that the majority of mismatches are due either to an inaccurate calculation of drainage area or
#imprecise or inaccurate coordinates in the ANEEL database.

#For large KWTOQ ratio (> 1000 electricty production (KW)/ annual average discharge (m3/s), 500 dams:
    #Make sure that dam wasn't snapped to smaller tributary
    #For KWTOQ > > 5000: Check that river name matches and that if in operation, could see reservoir or station.
    # Looked for comments (e.g. incorrect coordinates)
#For small KWTOQ ratio (< 10):
    #Make sure that small dams didn't get snapped to large rivers. did not find any anomaly. Could also be typos in production capacity.
#For dams 200-500 m (766 dams), manually inspect and snap either:
    #  already in operation and reasonable certainty that reservoir or diversion exists and/or original point sits on top of watercourse in imagery and/or brazilian hydrography
    #  original point sits on top of watercourse in imagery and/or brazilian hydrography
    # when in doubt, make sure that river name matches
#Check every dam > 30000 KW

#Check every dam that has DRENCOMP > 1.25 or DRENCOMP < 0.75
#Even when dam was topologically placed the same way in the network as on the national river network, we moved points
#to places where there was an obvious correspondence to the reported drainage area in the database (moved past or up a confluence, from tributary to mainstem, etc.).
#Many instances when drainage area did not match but the dam was topologically well places, there no stream in the
# surrounding area matching drainage area, often a reservoir or obvious waterfall sites was visible.

#Make a copy of damsnap_edit to avoid deleting by mistake
arcpy.CopyFeatures_management(damsnap_edit, os.path.join(datadir, damsnap_edit + 'copy'))

#Re-snap (to make sure that manual edits are properly intersecting network) and delete dams
arcpy.CopyFeatures_management(damsnap_edit, damsnap_editclean)
arcpy.Snap_edit(damsnap_editclean, [snap_env])
with arcpy.da.UpdateCursor(damsnap_editclean, ['manualsnap', 'DAMID']) as cursor:
    for row in cursor:
        if row[0]== '-1':
            print('Delete {}'.format(row[1]))
            cursor.deleteRow()

########################################################################################################################
# NETWORK PREPARATION
########################################################################################################################
#---- 6. Assign reaches to basin for each level, project them, and create geometric network ----
#Get geodesic length of full reaches prior to intersecting with basins
if not 'LENGTH_GEOFULL' in [f.name for f in arcpy.ListFields(snbr_clean)]:
    arcpy.AddGeometryAttributes_management(snbr_clean, 'LENGTH_GEODESIC', Length_Unit = "METERS")
    arcpy.AlterField_management(snbr_clean, 'LENGTH_GEO', new_field_name='LENGTH_GEOFULL', new_field_alias = 'LENGTH_GEOFULL')

#For each level
#for level in blevels:
level = '08'
print('Preparing network for level {}'.format(level))
#Intersect stream reaches with basins at multiple levels, leading to some reaches being divided in multiple segments
basinbrproj = os.path.join(dcigdb, "BasinATLAS_lev{}_v10_brproj".format(level))
outsn = os.path.join(dcigdb, 'streamATLAS_intersect_level{}'.format(level))

if not arcpy.Exists(outsn):
    print('Intersecting stream network with {}'.format(basinbrproj))
    arcpy.Intersect_analysis([snbr_clean, basinbrproj], out_feature_class=outsn, join_attributes='ALL')
else:
    print('{} already exists. Skipping...'.format(outsn))

arcpy.AddGeometryAttributes_management(outsn, 'LENGTH_GEODESIC', Length_Unit = 'METERS')

###Assign basin that reaches overlaps the most to all subsegments of reach
#Identify basin that each reach overlaps the most with and compute # of intersecting basins
basoverlap = defaultdict(list) #Output: [length of subsegment, HYBAS_ID, # of intersecting basins]
for row in arcpy.da.SearchCursor(outsn, ['REACH_ID', 'LENGTH_GEO', 'HYBAS_ID']):
    if row[0] in basoverlap:
        basoverlap[row[0]][2] += 1 #Add 1 to number of intersecting basins
        if row[1] > basoverlap[row[0]][0]:
            basoverlap[row[0]][0] = row[1]
            basoverlap[row[0]][1] = row[2]
    else:
        basoverlap[row[0]] = [row[1], row[2], 1]

#Assign it to the original stream network
arcpy.AddField_management(in_table = snbr_clean, field_name = 'HYBAS_ID{}'.format(level),
                          field_type=[f.type for f in arcpy.ListFields(outsn) if f.name == 'HYBAS_ID'][0])
with arcpy.da.UpdateCursor(snbr_clean, ['REACH_ID', 'LENGTH_GEOFULL', 'HYBAS_ID{}'.format(level)]) as cursor:
    for row in cursor:
        if row[0] in basoverlap:
            #Delete segments that only intersect one basin within brazil and are not fully contained within it
            #(deals with stubs that flow to basins outside of Brazil)
            if (basoverlap[row[0]][2] == 1 and round(basoverlap[row[0]][0],1) < round(row[1],1)):
                #row[2] = '-999.999'
                cursor.deleteRow()
            else:
                row[2] = basoverlap[row[0]][1]
                cursor.updateRow(row)

#----------------------------- Split network at dams -------------------------------------
print('Split stream segments where dams intersect them...')
arcpy.SplitLineAtPoint_management(snbr_clean, damsnap_editclean, outsnsplit, search_radius= '0.1 Meters')

#Create new unique id for each segment in network
outfield = 'SUBREACH_ID'
flist = [f.name for f in arcpy.ListFields(outsnsplit)]
if outfield not in flist:
    print('Create SUBREACH_ID field...')
    arcpy.AddField_management(outsnsplit, outfield, 'LONG')

edit.startEditing(False, True)
print('Compute SUBREACH_ID field...')
i=0
with arcpy.da.UpdateCursor(outsnsplit, [outfield]) as cursor:
    for row in cursor:
        row[0] = i
        cursor.updateRow(row)
        i += 1

edit.stopEditing(True)

#-------------- Edit dams that fall on the start node (upstream edge) of 1st order streams ------------------------------------------
#Find dangle points
print('Feature vertices to points...')
arcpy.FeatureVerticesToPoints_management(outsnsplit, '{}_vertices'.format(outsnsplit), point_location= 'DANGLE')

#Intersect
print('Intersect network vertices with dams...')
dangledam = os.path.join(dcigdb, 'damdangle_inters')
arcpy.Intersect_analysis(['{}_vertices'.format(outsnsplit), damsnap_editclean], dangledam,
                         join_attributes='ALL')

#Make list of dam-line pairs
damlinedict = dict([(row[0], row[1]) for row in arcpy.da.SearchCursor(dangledam, ['SUBREACH_ID', 'DAMID'])])

#Shift dam down line
min([row[0] for row in arcpy.da.SearchCursor(dangledam, ['LENGTH_GEOFULL'])]) #Check length of lines that dams sit on
#Get position of point intersection of 200-m buffer around point with line
dangledambuf = os.path.join(dcigdb, 'damdangle_buf')
arcpy.Buffer_analysis(dangledam, dangledambuf, '200 meters')
dangledambufline = os.path.join(dcigdb, 'damdangle_bufline')
arcpy.FeatureToLine_management(dangledambuf, dangledambufline, attributes='ATTRIBUTES')
dangledambufinters = os.path.join(dcigdb, 'damdangle_bufpoint')
arcpy.Intersect_analysis([outsnsplit, dangledambufline], dangledambufinters, output_type='POINT', join_attributes='ALL')
#Get XY coordinates of intersection of the buffer with the line that the dam sits on (rather than potentially other nearby line)
damlineshiftdict = dict([(row[2], row[3]) for row in
                         arcpy.da.SearchCursor(dangledambufinters, ['SUBREACH_ID', 'SUBREACH_ID_1', 'DAMID', 'SHAPE@XY'])
                         if row[0] == row[1]])
#Edit dam layer to shift dams
with arcpy.da.UpdateCursor(damsnap_editclean, ['DAMID', 'SHAPE@XY']) as cursor:
    for row in cursor:
        if row[0] in damlineshiftdict:
            row[1] = damlineshiftdict[row[0]]
            cursor.updateRow(row)

#----------------------------- Split network at dams again -------------------------------------
print('Split stream segments where dams intersect them...')
arcpy.SplitLineAtPoint_management(snbr_clean, damsnap_editclean, outsnsplit, search_radius= '0.010 Meters')

#Create new unique id for each segment in network
outfield = 'SUBREACH_ID'
flist = [f.name for f in arcpy.ListFields(outsnsplit)]
if outfield not in flist:
    print('Create SUBREACH_ID field...')
    arcpy.AddField_management(outsnsplit, outfield, 'LONG')

edit.startEditing(False, True)
print('Compute SUBREACH_ID field...')
i=0
with arcpy.da.UpdateCursor(outsnsplit, [outfield]) as cursor:
    for row in cursor:
        row[0] = i
        cursor.updateRow(row)
        i += 1
edit.stopEditing(True)

#Compute segment length
if 'LENGTH' not in flist:
    print('Compute segment length...')
    arcpy.AddGeometryAttributes_management(outsnsplit, 'LENGTH', 'meters')

#----------------------------- Create geometric network -------------------------------------
in_source_feature_classes = "{} SIMPLE_EDGE NO".format(os.path.split(outsnsplit)[1])
if not arcpy.Exists(geonetout):
    print('Create geometric network...')
    arcpy.CreateGeometricNetwork_management(netftdat, os.path.split(geonetout)[1], in_source_feature_classes)
else:
    print('{} already exists. Skipping...'.format(geonetout))

########################################################################################################################
# DCI ANALYSIS
########################################################################################################################
#Set flow direction of geometric network, a pre-requisite for tracing
arcpy.SetFlowDirection_management(geonetout, flow_option = "WITH_DIGITIZED_DIRECTION")

#For each segment, identify connected reaches with dams as barriers
#Get start point of each reach
arcpy.FeatureVerticesToPoints_management(outsnsplit, netp, point_location='START')
segdic = defaultdict(list) #dictionary to contain segment ID as key and SUBREACH_ID of corresponding reaches as values
#Iterate over each reach
print('Find connected reaches...')
#################### FUTURE EDIT TO DO: WRITE DIRECTLY DIC WITH KEYS AS SUBREACH_id AND VALUE AS SEGID TO MAKE IT FASTER - FINE FOR NOW
with arcpy.da.SearchCursor(netp, ['SUBREACH_ID', 'SHAPE@']) as cursor:
    i=0
    for row in cursor:
        if row[0] not in set(itertools.chain.from_iterable(segdic.values())): #If reach not already part of a segment
            print(i)
            #Select all connected reaches using dams a barriers (this includes stopping reaches where dam is though)
            arcpy.TraceGeometricNetwork_management(in_geometric_network= geonetout,
                                                   out_network_layer = "connectlyr",
                                                   in_flags = row[1], in_barriers= damsnap_editclean,
                                                   in_trace_task_type= "FIND_CONNECTED",
                                                   in_trace_ends = 'NO_TRACE_ENDS')
            connect_tr = "connectlyr/{}".format(os.path.split(outsnsplit)[1])

            #Add reaches with barriers on them that are connected
            arcpy.TraceGeometricNetwork_management(in_geometric_network= geonetout,
                                                   out_network_layer = "connectstoplyr",
                                                   in_flags = row[1], in_barriers=damsnap_editclean,
                                                   in_trace_task_type= "FIND_CONNECTED",
                                                   in_trace_ends = 'TRACE_ENDS')
            connectstop_tr = "connectstoplyr/{}".format(os.path.split(outsnsplit)[1])

            #Populate dictionary with SUBREACH_ID for that segment ID
            with arcpy.da.SearchCursor(connect_tr, ["SUBREACH_ID"]) as cursor_connect:
                for row_connect in cursor_connect:
                    segdic[i].append(row_connect[0])
            with arcpy.da.SearchCursor(connectstop_tr, ["SUBREACH_ID"]) as cursor_connectstop:
                for row_connectstop in cursor_connectstop:
                    segdic[i].append(row_connectstop[0])

            #Check if the segment is empty (sometimes happens for single reaches between two dams) and fill it in with the corresponding reach
            if len(segdic[i]) == 0:
                countout = arcpy.GetCount_management("connectstoplyr/{}_Junctions".format(os.path.split(geonetout)[1]))
                if int(countout.getOutput(0)) == 1:
                    print('Assign SEGID {0} to reach {1}'.format(i, row[0]))
                    segdic[i].append(row[0])
            i+=1

#Reverse dictionary to assign a segment ID to each reach
invsegdic = {v: k for k, vall in segdic.iteritems() for v in vall}

# #### In case one wants to save dictionary to file
import cPickle as pickle
outpickle = os.path.join(resdir, 'invsegdic')
pickle.dump(invsegdic, open(outpickle, "wb" ) )
invsegdic = pickle.load(open(outpickle, "rb" ) )

#Assign segment ID to each reach in network
print('Assign segment ID to reach reach...')
if 'SEGID' not in [f.name for f in arcpy.ListFields(outsnsplit)]:
    arcpy.AddField_management(outsnsplit, 'SEGID', 'LONG')
edit.startEditing(False, True)
with arcpy.da.UpdateCursor(outsnsplit, ['SUBREACH_ID', 'SEGID']) as cursor:
    for row in cursor:
        try:
            row[1] = invsegdic[row[0]]
            cursor.updateRow(row)
        except:
            print('{} not in list of reaches'.format(row[0]))
edit.stopEditing(True)

#Dissolve based on SEGID and HYBAS_ID08 (to get at higher or lower level of basins HYBAS_ID06, or HYBAS_ID10, dissolve by these fields)
print('Dissolve network based on segment ID and basin ID...')
arcpy.Dissolve_management(outsnsplit, net_diss, dissolve_field = ['HYBAS_ID08', 'SEGID'])

#Link dams with upstream and downstream reach thanks to last and first point of SHAPE@ (use int(shape.X) and Y to deal with points sometimes being a few cm off)
print('Identify upstream and downstream segment for every dam...')
arcpy.FeatureVerticesToPoints_management(net_diss, netdanglep, point_location='DANGLE')
arcpy.FeatureVerticesToPoints_management(net_diss, netendp, point_location='END')

#Dictionary containing rounded coordinates and corresponding SEGID for start vertex of every segment (to identify downstream segments)
########
#Future edit: Exclude dangle points to avoid errors when dams are too close to the start vertex of first order stream (e.g. D1430 Ipueiras)
#However, issues of resolution/tolerance made a point lower down a dangle in this case. Not sure how to solve. Correct manually in R code
#danglepset = {(round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1)) for row in arcpy.da.SearchCursor(netdanglep, ['SHAPE@'])}
####
firstpdict = defaultdict(list)
for row in arcpy.da.SearchCursor(netp, ['SHAPE@', 'SUBREACH_ID']):
    outcoor = (round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))
    #if outcoor not in danglepset:
    firstpdict[outcoor].append(invsegdic[row[1]])

#Dictionary containing rounded coordinates and corresponding SEGID for end vertex of every segment (to identify upstream segments)
lastpdict = defaultdict(list)
for row in arcpy.da.SearchCursor(netendp, ['SHAPE@', 'SEGID']):
    lastpdict[(round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))].append(row[1])

#Create UpSeg and DownSeg fields
damfieldlist = [f.name for f in arcpy.ListFields(damsnap_editclean)]
if 'UpSeg' not in damfieldlist:
    arcpy.AddField_management(damsnap_editclean, 'UpSeg', 'LONG')
if 'downSeg' not in damfieldlist:
    arcpy.AddField_management(damsnap_editclean, 'DownSeg', 'LONG')

#Assign UpSeg and DownSeg for each dam by overlapping it with start and end vertices
#Manual corrections that can be improved upon automatically in later edits
bugdictstart = {}

with arcpy.da.UpdateCursor(damsnap_editclean, ['SHAPE@', 'DownSeg', 'UpSeg', 'DAMID']) as cursor:
    for row in cursor:
        rowp = (round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))
        if rowp in firstpdict.keys():
            row[1] = firstpdict[rowp][0]
        if row[3] in bugdictstart: #If known glitch
            row[1] = bugdictstart[row[3]]
        if rowp in lastpdict.keys():
            for p in lastpdict[rowp]: #Deal with bug when downsegment has an end point overlapping with that of the upsegment (not sure why)
                if p != firstpdict[rowp][0]:
                    row[2] = p
        cursor.updateRow(row)

#Spatial join dams to network to get SUBREACH_ID
print('Spatial join dam to network datasets')
arcpy.SpatialJoin_analysis(damsnap_editclean, outsnsplit, damsnap_editcleanjoin, join_operation='JOIN_ONE_TO_ONE',
                           match_option = 'INTERSECT', search_radius = '0.001 meters')

#Recreate unique DAMID as manual edits created some duplicates (shouldn't occur when re-running code/edits)
with arcpy.da.UpdateCursor(damsnap_editcleanjoin, [arcpy.Describe(damsnap_editcleanjoin).OIDFieldName, 'DAMID']) as cursor:
    for row in cursor:
        row[1] = 'D{}'.format(row[0])
        cursor.updateRow(row)

#Add X and Y to dams (in WGS84)
arcpy.AddGeometryAttributes_management(Input_Features=damsnap_editcleanjoin, Geometry_Properties=['POINT_X_Y_Z_M'],
                                       Coordinate_System=arcpy.SpatialReference(4326))

#Export dam and network attributes to tables
arcpy.CopyRows_management(damsnap_editcleanjoin, dam_tab)
arcpy.CopyRows_management(net_diss, netseg_tab)