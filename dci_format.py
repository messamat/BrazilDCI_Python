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
rootdir = 'F:/Mathis/brazilDCI'
blevels = ['04', '06', '08', '10']
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')
crsSIRGAS = arcpy.SpatialReference(5880) #SIRGAS_2000_Brazil_Polyconic (ESPG: 5880)
streamatlas = os.path.join(datadir, 'HydroATLAS/StreamATLAS_v01.gdb/StreamATLAS_v01')
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

basinbr4 = os.path.join(dcigdb, "BasinATLAS_lev04_v01_br")
netproj4 = os.path.join(netftdat, 'streamATLAS_proj_level04')
netproj10 = os.path.join(netftdat, 'streamATLAS_proj_level10')
geonet4 = os.path.join(netftdat, 'streamATLAS_geonet_level4')
geonet10 = os.path.join(netftdat, 'streamATLAS_geonet_level10')

dams_original = os.path.join(datadir, 'DamsBrazil/Aproveitamento_Hidrelétricos_AHE.shp')
dams_originalgdb = os.path.join(dcigdb, 'Aproveitamento_Hidreletricos_AHE')
damsbuf300 = os.path.join(dcigdb, 'damsbuf300')
damsbufinter = os.path.join(dcigdb, 'damsbufinter')
damsel = os.path.join(dcigdb, 'damfull_simplified_300m')
damsnap = os.path.join(dcigdb, 'damfull_snapped')
damsnap_edit = os.path.join(dcigdb, 'damfull_snappededit')
damsnap_join = os.path.join(dcigdb, 'damfull_join')

outsnsplit = os.path.join(netftdat, 'streamATLAS_split')
geonetout = os.path.join(netftdat,'streamATLAS_geonet')

netp = os.path.join(dcigdb, 'netpointstart')
netendp = os.path.join(dcigdb, 'netpointend')
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
    #---- 1. Filter basins and rivers of Brazil ----
    basin = os.path.join(datadir, "HydroATLAS/BasinATLAS_v01.gdb/BasinATLAS_lev{}_v01".format(level))
    basinbr = os.path.join(dcigdb, "{}_br".format(os.path.split(basin)[1]))
    basinbrinter = os.path.join(dcigdb, "{}_brinter".format(os.path.split(basin)[1]))
    if not arcpy.Exists(basinbr):
        print('Filter BasinATLAS for Brazil for level {}...'.format(level))
        if not 'AREA_GEOFULL' in [f.name for f in arcpy.ListFields(basin)]:
            arcpy.AddGeometryAttributes_management(basin, 'AREA_GEODESIC', Area_Unit = "SQUARE_KILOMETERS")
            arcpy.AlterField_management(basin, 'AREA_GEO', new_field_name='AREA_GEOFULL', new_field_alias = 'AREA_GEOFULL')
        arcpy.Intersect_analysis([basin, br_bound], basinbrinter)
        arcpy.AddGeometryAttributes_management(basinbrinter, 'AREA_GEODESIC', Area_Unit = 'SQUARE_KILOMETERS')
        arcpy.MakeFeatureLayer_management(basinbrinter, 'basinbrinterlyr',
                                          where_clause= '(AREA_GEO/AREA_GEOFULL) > 0.01')
        basinbr4list = [row[0] for row in arcpy.da.SearchCursor('basinbrinterlyr', ['HYBAS_ID'])]
        arcpy.MakeFeatureLayer_management(basin, 'basinbr4lyr',
                                          where_clause= 'HYBAS_ID IN {}'.format(tuple(basinbr4list)))
        arcpy.CopyFeatures_management('basinbr4lyr', basinbr)

    #---- 4. Create a shapefile removing coastal drainages and project basins----
    basinnocoast =  os.path.join(dcigdb, "{}_NonCoastal".format(basinbr))
    basinnocoastproj = os.path.join(dcigdb, "{}proj".format(basinnocoast))
    if not arcpy.Exists(basinnocoastproj):
        print('Remove coastal drainages from BasinATLAS for level {}...'.format(level))
        arcpy.MakeFeatureLayer_management(basinbr, 'basinnocoastlyr')
        arcpy.SelectLayerByLocation_management('basinnocoastlyr', 'INTERSECT', ngacoast,
                                               invert_spatial_relationship='INVERT')
        arcpy.CopyFeatures_management('basinnocoastlyr', basinnocoast)
        arcpy.Project_management(basinnocoast, basinnocoastproj, crsSIRGAS)

#---- 2. Subset StreamATLAS based on level 4 basins that intersect Brazil ----
arcpy.Clip_analysis(streamatlas, basinbr4, snbr_original)

#---- 3. Create a shapefile only for dams in operation ----


#-----5. Clean and project river network + project dams ----
#Identify and remove the bugs in the streamATLAS network
#Project dam dataset and import it into gdb
if not arcpy.Exists(dams_originalgdb):
    arcpy.Project_management(dams_original, dams_originalgdb, crsSIRGAS)

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
with arcpy.da.UpdateCursor(dams_originalgdb, [arcpy.Describe(damsel).OIDFieldName, 'DAMID']) as cursor:
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
snap_env = [snbr_clean, "EDGE", "500 Meters"]
arcpy.Snap_edit(damsnap, [snap_env])

#Check those dams > 500 m from river network with imagery. Either move them to the network or delete them.
#Create field to know whether dam was checked (0) or edited (1), otherwise 'null'
arcpy.CopyFeatures_management(damsnap, damsnap_edit)
arcpy.AddField_management(damsnap_edit, 'manualsnap', 'TEXT')
damdel = {} #list of dams to delete
#Delete dams
with arcpy.da.UpdateCursor(damsnap_edit, [arcpy.Describe(damsel).OIDFieldName]) as cursor:
    for row in cursor:
        if row[0] in damdel:
            cursor.deleteRow()

#For now, delete all those < 500 m, refine later
arcpy.MakeFeatureLayer_management(damsnap, 'damsnaplyr', where_clause='NEAR_DIST < 500')
arcpy.CopyFeatures_management('damsnaplyr', damsnap_edit)

########################################################################################################################
# NETWORK PREPARATION
########################################################################################################################
#---- 6. Assign reaches to basin for each level, project them, and create geometric network ----
#Create feature dataset
#pr.XYTolerance = 0.01
if not arcpy.Exists(netftdat):
    arcpy.CreateFeatureDataset_management(dcigdb, os.path.split(netftdat)[1], crsSIRGAS) #pr
#arcpy.Describe(netftdat).SpatialReference.XYTolerance

#For each level
for level in blevels:
    #Intersect stream reaches with basins at multiple levels, leading to some reaches being divided in multiple segments
    basinnocoast = os.path.join(dcigdb, "BasinATLAS_lev{}_v01_br_NonCoastalproj".format(level))
    outsn = os.path.join(dcigdb, 'streamATLAS_intersect_level{}'.format(level))
    if not arcpy.Exists(outsn):
        print('Intersecting stream network with {}'.format(basinnocoast))
        arcpy.Intersect_analysis([snbr_clean, basinnocoast], out_feature_class= outsn, join_attributes='ALL')
    else:
        print('{} already exists. Skipping...'.format(outsn))

    ###Assign basin that reaches overlaps the most to all subsegments of reach
    #Identify basin that each reach overlaps the most with
    basoverlap = defaultdict(list)
    for row in arcpy.da.SearchCursor(outsn, ['REACH_ID', 'Shape_Length', 'HYBAS_ID']):
        if row[0] in basoverlap:
            if row[1] > basoverlap[row[0]][0]:
                basoverlap[row[0]] = [row[1], row[2]]
        else:
            basoverlap[row[0]] = [row[1], row[2]]
    #Assign it to the original stream network
    arcpy.AddField_management(in_table = snbr_clean, field_name = 'HYBAS_ID{}'.format(level),
                              field_type=[f.type for f in arcpy.ListFields(outsn) if f.name == 'HYBAS_ID'][0])
    with arcpy.da.UpdateCursor(snbr_clean, ['REACH_ID', 'HYBAS_ID{}'.format(level)]) as cursor:
        for row in cursor:
            if row[0] in basoverlap:
                row[1] = basoverlap[row[0]][1]
                cursor.updateRow(row)

#----------------------------- Split network at dams -------------------------------------
if not arcpy.Exists(outsnsplit):
    print('Split stream segments where dams intersect them...')
    arcpy.SplitLineAtPoint_management(snbr_clean, damsnap_edit, outsnsplit, search_radius= '0.1 Meters')

#-------------- Edit dams that fall on the start node (upstream edge) of 1st order streams ------------------------------------------
#Find dangle points
arcpy.FeatureVerticesToPoints_management(outsnsplit, '{}_vertices'.format(outsnsplit), point_location= 'DANGLE')

#Intersect
dangledam = os.path.join(dcigdb, 'damdangle_inters')
arcpy.Intersect_analysis(['{}_vertices'.format(outsnsplit), damsnap_edit], dangledam,
                         join_attributes='ALL')

#Make list of dam-line pairs
damlinedict = dict([(row[0], row[1]) for row in arcpy.da.SearchCursor(dangledam, ['SUBREACH_ID', 'DAMID'])])

#Shift dam down line
min([row[0] for row in arcpy.da.SearchCursor(dangledam, ['LENGTH'])]) #Check length of lines that dams sit on
#Get position of point intersection of 50-m buffer around point with line
dangledambuf = os.path.join(dcigdb, 'damdangle_buf')
arcpy.Buffer_analysis(dangledam, dangledambuf, '50 meters')
dangledambufline = os.path.join(dcigdb, 'damdangle_bufline')
arcpy.FeatureToLine_management(dangledambuf, dangledambufline, attributes='ATTRIBUTES')
dangledambufinters = os.path.join(dcigdb, 'damdangle_bufpoint')
arcpy.Intersect_analysis([outsnsplit, dangledambufline], dangledambufinters, output_type='POINT', join_attributes='ALL')
#Get XY coordinates of intersection of the buffer with the line that the dam sits on (rather than potentially other nearby line)
damlineshiftdict = dict([(row[2], row[3]) for row in
                         arcpy.da.SearchCursor(dangledambufinters, ['SUBREACH_ID', 'SUBREACH_ID_1', 'DAMID', 'SHAPE@XY'])
                         if row[0] == row[1]])
#Edit dam layer to shift dams
with arcpy.da.UpdateCursor(damsnap_edit, ['DAMID', 'SHAPE@XY']) as cursor:
    for row in cursor:
        if row[0] in damlineshiftdict:
            row[1] = damlineshiftdict[row[0]]
            cursor.updateRow(row)

#----------------------------- Split network at dams again -------------------------------------
if not arcpy.Exists(outsnsplit):
    print('Split stream segments where dams intersect them...')
    arcpy.SplitLineAtPoint_management(snbr_clean, damsnap_edit, outsnsplit, search_radius= '0.001 Meters')

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
    "DElete first"

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
#################### EDIT TO DO: WRITE DIRECTLY DIC WITH KEYS AS SUBREACH_id AND VALUE AS SEGID
with arcpy.da.SearchCursor(netp, ['SUBREACH_ID', 'SHAPE@']) as cursor:
    i=0
    for row in cursor:
        if row[0] not in set(itertools.chain.from_iterable(segdic.values())): #If reach not already part of a segment
            print(i)
            #Select all connected reaches using dams a barriers (this includes stopping reaches where dam is though)
            arcpy.TraceGeometricNetwork_management(in_geometric_network= geonetout,
                                                   out_network_layer = "connectlyr",
                                                   in_flags = row[1], in_barriers= damsnap_edit,
                                                   in_trace_task_type= "FIND_CONNECTED",
                                                   in_trace_ends = 'NO_TRACE_ENDS')
            connect_tr = "connectlyr/{}".format(os.path.split(outsnsplit)[1])

            #Add reaches with barriers on them that are connected
            arcpy.TraceGeometricNetwork_management(in_geometric_network= geonetout,
                                                   out_network_layer = "connectstoplyr",
                                                   in_flags = row[1], in_barriers=damsnap_edit,
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

# #In case one wants to save dictionary to file
# import cPickle as pickle
# outpickle = os.path.join(resdir, 'invsegdic')
# pickle.dump(invsegdic, open(outpickle, "wb" ) )
# invsegdic = pickle.load(open(outpickle, "rb" ) )

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

#Dissolve based on SEGID and HYBAS_ID08
print('Dissolve network based on segment ID and basin ID...')
arcpy.Dissolve_management(outsnsplit, net_diss, dissolve_field = ['HYBAS_ID08', 'SEGID'])

#Link dams with upstream and downstream reach thanks to last and first point of SHAPE@ (use int(shape.X) and Y to deal with points sometimes being a few cm off)
print('Identify upstream and downstream segment for every dam...')
arcpy.FeatureVerticesToPoints_management(net_diss, netendp, point_location='END')
lastpdict = defaultdict(list)
for row in arcpy.da.SearchCursor(netendp, ['SHAPE@', 'SEGID']):
    lastpdict[(round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))].append(row[1])

firstpdict = defaultdict(list)
for row in arcpy.da.SearchCursor(netp, ['SHAPE@', 'SUBREACH_ID']):
    firstpdict[(round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))].append(invsegdic[row[1]])

damfieldlist = [f.name for f in arcpy.ListFields(damsnap_edit)]
if 'UpSeg' not in damfieldlist:
    arcpy.AddField_management(damsnap_edit, 'UpSeg', 'LONG')
if 'downSeg' not in damfieldlist:
    arcpy.AddField_management(damsnap_edit, 'DownSeg', 'LONG')

with arcpy.da.UpdateCursor(damsnap_edit, ['SHAPE@', 'DownSeg', 'UpSeg', 'DAMID']) as cursor:
    for row in cursor:
        rowp = (round(row[0].centroid.X, 1), round(row[0].centroid.Y, 1))
        if rowp in firstpdict.keys():
            row[1] = firstpdict[rowp][0]
        if rowp in lastpdict.keys():
            for p in lastpdict[rowp]: #Deal with bug when downsegment has an end point overlapping with that of the upsegment (not sure why)
                if p != firstpdict[rowp][0]:
                    row[2] = p
        cursor.updateRow(row)

#Get HYBASID_08
#Spatial join dams to network to get SUBREACH_ID
arcpy.SpatialJoin_analysis(damsnap_edit, outsnsplit, damsnap_join, join_operation='JOIN_ONE_TO_ONE',
                           match_option = 'INTERSECT', search_radius = '0.001 meters')

#Export dam and network attributes to tables
arcpy.CopyRows_management(damsnap_join, dam_tab)
arcpy.CopyRows_management(net_diss, netseg_tab)

########################################################################################################################
# EXTRA STUFF BY MATHIS
########################################################################################################################
# #Dissolve network by HYBAS_ID
# arcpy.env.workspace = dcigdb
# [f.name for f in arcpy.ListFields(netsub)]
# arcpy.Dissolve_management(netsub, 'netsubdiss', dissolve_field='PFAF_ID', statistics_fields=[['LENGTH_GEO', 'SUM']],
#                           multi_part='MULTI_PART', unsplit_lines='DISSOLVE_LINES')
#

#### Delete upstream stream segments that dangle in adjacent basin after intersecting
    # splitdic = defaultdict(int)
    # reachnumbas = defaultdict(int)
    # for row in arcpy.da.SearchCursor(outsn, ['REACH_ID', 'HYBAS_ID']):
    #     splitdic[row[0]] += 1 #Identify reaches that were split when intersecting with basins (as > 1 feature with same REACH_ID)
    #     reachnumbas[row[1]] += 1 #Count number of reaches per basin
    #
    # arcpy.FeatureVerticesToPoints_management(outsn, '{}_vertices'.format(outsn), point_location= 'DANGLE')
    # #Get ID of all streams with dangle points if at least one dangle point is > 0.36 km from mouth (half diagonal of a pixel)
    # #and if more than one line in basin (if not, it's just that this is the only stream segment which starts and end basin)
    # dangledict = {}
    # for row in arcpy.da.SearchCursor('{}_vertices'.format(outsn), ['REACH_ID', 'ORIG_FID', 'rld_km_pva', 'HYBAS_ID']):
    #     if row[2] > 0.36 and reachnumbas[row[3]]>1:
    #         dangledict[row[0]] = row[1]
    #
    # dellist = []
    # with arcpy.da.UpdateCursor(outsn, ['REACH_ID', arcpy.Describe(outsn).OIDFieldName, 'Shape_Length']) as cursor:
    #     for row in cursor:
    #         if splitdic[row[0]] > 1: #If reach was split
    #             if (row[0] in dangledict and row[1] == dangledict[row[0]]) or row[2] < 0.05: #And reach segment has a danglepoint or is less than 50 meters (to get those that cross corners)
    #                 dellist.append(row[1])
    #
    # #dellist = [dangledict[k] for k,v in splitdic.iteritems() if v == 2 and k in dangledict]
    # arcpy.MakeFeatureLayer_management(outsn, 'updangle',
    #                                   '{0} IN {1}'.format(arcpy.Describe(outsn).OIDFieldName, str(tuple(dellist))))
    # arcpy.CopyFeatures_management('updangle', os.path.join(dcigdb, 'check'))
    #
    #
    # for row in arcpy.da.UpdateCursor(outsn, ['REACH_ID'])
    #
    #
    # duplilist = [id for id in [row[0] for row in arcpy.da.SearchCursor(outsn, ['REACH_ID'])]]
    #

########################################################################################################################
# EXTRA STUFF FOR BAT ANALYSIS
########################################################################################################################

################# NETWORK PREPARATION #####################
#---- 10. Export data of “streamATLAS_ proj_LevelXX.shp” to a layer to be modified by BAT
#a.	“BAT_network_Level08.shp”
#b.	“BAT_network_Level06.shp”
#c.	“BAT_network_Level04.shp”
#d.	“BAT_network_Level10.shp”

#----11) Create a column with basin IDs to run in BAT (RegionX). It should be a string (text field).
#  Use the field calculator to fill out the column with the HYBAS_ID.
# a.	“BAT_network_Level08.shp”
# b.	“BAT_network_Level06.shp”
# c.	“BAT_network_Level04.shp”
# d.	“BAT_network_Level10.shp”

#----12)	Optimize the analysis only running it for basins with dams.
# First, select manually the basins that have multiple networks (Coastal drainages) in the Level 4.
# Save as the layer as “BasinATLAS_level04_Brazil_NonCoastal.shp”, which will be used for removing the coastal
# drainages from all the other levels (keep the same basins for comparison). Second, select the basins with no
# coastal drainages using the “Select by Location” and the raw BasinATLAS for all levels as the target.
# se the source: “BasinATLAS_level04_Brazil_NonCoastal”. Export data and save all these resulting layers as
#  “BasinATLAS_levelXX_Brazil_NonCoastal.shp”.
# Third, use the “Select by location” on the “BasinATLAS_levelXX_Brazil_NonCoastal” layer and the source
#  layer “DamAll_Proposed.shp”. Export data and save as “BasinsLX_WithDams.shp”.
# Finally, clip the “BAT_network_LevelXX.shp” layer using the “BasinsLX_WithDams.shp” as the clipping feature.
#a.	“BAT_networkX_WithDams.shp”

# 13)	Export the data “BAT_networkX_WithDams.shp” for future edition by the BAT analysis
# a.	“BAT_RUN_Net10.shp”
# b.	“BAT_RUN_Net08.shp”
# c.	“BAT_RUN_Net06.shp”
# d.	“BAT_RUN_Net04.shp”

################ DAM DATASET PREPARATION #####################
# 17)	Remove coastal drainages from “DamFull_Simplified_debug_FINAL.shp”.
# “Geoprocessing” -> “Clip”. Based on level04 ATLAS with non-coastal basins.
# a.	“DammFull_Simplified_NonCoastal.shp”
