__author__ = 'Mathis Messager'

import arcpy
import os
from collections import defaultdict


#Set analysis parameters
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = True

#Input data and variables
rootdir = 'F:/brazilDCI'
blevels = ['04', '06', '08', '10']
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')
crsSIRGAS = arcpy.SpatialReference(4674)

#Output variables
snbr_original = os.path.join(resdir, 'streamATLAS_BrazilOriginal.shp')
dcigdb = os.path.join(resdir, 'dci.gdb')
if not os.path.exists(dcigdb):
    arcpy.CreateFileGDB_management(resdir, 'dci')
snbr_debug = os.path.join(dcigdb, 'streamATLAS_Brazil_debug201904')
netftdat = os.path.join(dcigdb, 'streamATLAS_geonetwork')
netproj4 = os.path.join(netftdat, 'streamATLAS_proj_level04')
netproj10 = os.path.join(netftdat, 'streamATLAS_proj_level10')

dams_original = os.path.join(datadir, 'DamsBrazil/Aproveitamento_Hidrelétricos_AHE.shp')
dams_originalgdb = os.path.join(dcigdb, 'Aproveitamento_Hidreletricos_AHE')
damsbuf300 = os.path.join(dcigdb, 'damsbuf300')
damsbufinter = os.path.join(dcigdb, 'damsbufinter')
damsel = os.path.join(dcigdb, 'damfull_simplified_300m')
damsnap = os.path.join(dcigdb, 'damfull_snapped')
damsnap_edit = os.path.join(dcigdb, 'damfull_snappededit')

########################################################################################################################
# GENERAL DATA PRE-FORMATTING
########################################################################################################################
#---- 1. Filter basins and rivers of Brazil ----

#---- 2. Create a shapefile only for dams in operation ----

#---- 3. Snap the dams to river segments ----

#---- 4. Create a shapefile removing coastal drainages ----

########################################################################################################################
# BARRIER PREPARATION
########################################################################################################################
#14) Clean the dam dataset
# (“Aproveitamento_Hidrelétricos_AHE.shp” – small and large/current and future).
# Merge centrals (e.g. Santa Lucia I and Santa Lucia II) and remove old ones (e.g. Dallasta, São Domingos and Ludesa).

#Import dataset into gdb
arcpy.CopyFeatures_management(dams_original, dams_originalgdb)

#Create unique dam ID
arcpy.AddField_management(dams_originalgdb, 'DAMID', 'TEXT')
with arcpy.da.UpdateCursor(dams_originalgdb, [arcpy.Describe(damsel).OIDFieldName, 'DAMID']) as cursor:
    for row in cursor:
        row[1] = 'D{}'.format(row[0])
        cursor.updateRow(row)

#create a buffer layer (300m) of the dams
buffdist = 300
arcpy.Buffer_analysis(dams_originalgdb, damsbuf300, '{} meters'.format(buffdist), line_side='FULL',
                      dissolve_option='NONE', method='GEODESIC')

# Second, identify what dam are within 300 m of each other and only keep the larger ones
arcpy.Intersect_analysis([dams_originalgdb, damsbuf300], damsbufinter, join_attributes='ALL')
arcpy.Sort_management(damsbufinter, damsbufinter + 'sort', sort_field=[['POT_KW', 'DESCENDING']])

[f.name for f in arcpy.ListFields(damsbufinter + 'sort')]
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

#Check that dataset is in SIRGAS2000 coord system
if arcpy.Describe(damsel).SpatialReference.name == crsSIRGAS.name:
    print('{} is already in SIRGAS 2000, no need to project...'.format(damsel))

# 16) Snap the dams with the streamATLAS-debug networks
#Check distance from dams to river network
arcpy.CopyFeatures_management(damsel, damsnap)
arcpy.Near_analysis(damsnap, netproj4, location = 'LOCATION', angle = 'NO_ANGLE', method = 'GEODESIC')

#Snap
snap_env = [netproj4, "EDGE", "500 Meters"]
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
#---- 5.	Identify and remove the bugs in the streamATLAS network (“streamATLAS_BrazilOriginal.shp” keeps untouched ----
# as a reference). They were identified during previous tests using BAT. - CREATE TOPOLOGICAL RULES TO AUTOMATICALLY DETECT THAT
outids = (60339808, 60340477, 60492207, 60571388, 60930843, 60955187, 60958534)
[[f.name, f.type] for f in arcpy.ListFields(snbr_original)]
arcpy.MakeFeatureLayer_management(snbr_original, 'sn_originallyr', "REACH_ID NOT IN {}".format(str(outids)))
arcpy.CopyFeatures_management('sn_originallyr', snbr_debug)

#---- 6. Divide the reaches that cross more than one basin, project them, and create geometric network ----
#Create feature dataset
#pr.XYTolerance = 0.01
if not arcpy.Exists(netftdat):
    arcpy.CreateFeatureDataset_management(dcigdb, os.path.split(netftdat)[1], crsSIRGAS) #pr
#arcpy.Describe(netftdat).SpatialReference.XYTolerance

#For each level
for level in blevels:
    #Intersect stream reaches with basins at multiple levels, leading to some reaches being divided in multiple segments
    inatlas = os.path.join(resdir, 'BasinATLAS_level{}_Brazil.shp'.format(level))
    outsn = os.path.join(dcigdb, 'streamATLAS_intersect_level{}'.format(level))
    if not arcpy.Exists(outsn):
        print('Intersecting stream network with {}'.format(inatlas))
        arcpy.Intersect_analysis([snbr_debug, inatlas], out_feature_class= outsn, join_attributes='ALL')
    else:
        print('{} already exists. Skipping...'.format(outsn))

    #Split network at dams
    outsnsplit = os.path.join(dcigdb, 'streamATLAS_split_level{}'.format(level))
    if not arcpy.Exists(outsnsplit):
        print('Split stream segments where dams intersect them...')
        arcpy.SplitLineAtPoint_management(outsn, damsnap_edit, outsnsplit, search_radius= '0.1 Meters')

    #Create new unique id for each segment in network
    outfield = 'SEGID'
    flist = [f.name for f in arcpy.ListFields(outsnsplit)]
    if outfield not in flist:
        print('Create SEGID field...')
        arcpy.AddField_management(outsnsplit, outfield, 'LONG')

    print('Compute SEGID field...')
    i=0
    with arcpy.da.UpdateCursor(outsnsplit, [outfield]) as cursor:
        for row in cursor:
            row[0] = i
            cursor.updateRow(row)
            i += 1

    #Project stream network
    print('Project to {}...'.format(crsSIRGAS.name))
    outproj = os.path.join(netftdat, 'streamATLAS_proj_level{}'.format(level))
    if not arcpy.Exists(outproj):
        arcpy.Project_management(outsnsplit, outproj, crsSIRGAS)
    else:
        print('{} already exists. Skipping...'.format(outproj))

    #Compute segment length
    if 'LENGTH_GEO' not in flist:
        print('Compute segment length...')
        arcpy.AddGeometryAttributes_management(outproj, 'LENGTH_GEODESIC', 'meters')

    #Create geometric network
    in_source_feature_classes = "{} SIMPLE_EDGE NO".format(os.path.split(outproj)[1])
    geonetout = 'streamATLAS_geonet_level{}'.format(level)
    if not arcpy.Exists(os.path.join(netftdat, geonetout)):
        print('Create geometric network...')
        arcpy.CreateGeometricNetwork_management(netftdat, geonetout, in_source_feature_classes)
    else:
        print('{} already exists. Skipping...'.format(geonetout))

##############################################################################################################
#ENSURE TOPOLOGICAl INTEGRITY

#Flag dams that fall on the upstream edge of 1st order streams
#Locate dams along geometric network reaches
routesub = os.path.join(dcigdb, "net10routes")
arcpy.CreateRoutes_lr(netproj10, "SEGID", routesub, "LENGTH")

postab = os.path.join(resdir, "damssegpos.dbf")
arcpy.LocateFeaturesAlongRoutes_lr(in_features=damsnap_edit, in_routes=routesub, route_id_field="SEGID",
                                   radius_or_tolerance='0.1 Meters', out_table=postab,
                                   out_event_properties= "SEGID POINT up_dist", in_fields = "FIELDS")

#Flag segments that split


########################################################################################################################
# DCI ANALYSIS
########################################################################################################################
#Start with subset of netproj 10
netsub = os.path.join(netftdat, 'netproj10_sub')
arcpy.AddField_management(netproj10, 'PFAF_ID4', 'LONG')

edit = arcpy.da.Editor(dcigdb)
edit.startEditing(False, True)

with arcpy.da.UpdateCursor(netproj10, ['PFAF_ID4', 'PFAF_ID']) as cursor:
    for row in cursor:
        row[0] = int(str(row[1])[0:4])
        cursor.updateRow(row)
arcpy.MakeFeatureLayer_management(netproj10, 'netproj10_lyr', where_clause='PFAF_ID4 = 6428')
#arcpy.GetCount_management('netproj10_lyr')
arcpy.CopyFeatures_management('netproj10_lyr', netsub)
edit.stopEditing(True)

arcpy.CreateGeometricNetwork_management(netftdat, 'netsubgeo',
                                        in_source_feature_classes = "{} SIMPLE_EDGE NO".format(os.path.split(netsub)[1]))
netsubgeo = os.path.join(netftdat, 'netsubgeo')


##################### IDENTIFY ALL UPSTREAM DOWN###############
#Create a list of ComID that dams sit on
dams_comid = []
with arcpy.da.SearchCursor(dams_rout, ["FLComID"]) as cursor_d:
    for row_d in cursor_d:
        dams_comid.append(row_d[0])

#Set flow direction of geometric network, a pre-requisite for tracing
arcpy.SetFlowDirection_management(network, flow_option = "WITH_DIGITIZED_DIRECTION")

#Iterate through gages
gages_sel = "gagetrim50k_o10y_netjoin"

arcpy.MakeFeatureLayer_management(gages_sel, "gage_sel_lyr")
gage_up_dic = defaultdict(list)
fields = arcpy.ListFields(gages_sel)
for f in fields:
    print f.name
#Iterate through gages (as the "TraceGeometricNetwork" tool does not allow for a tolerance parameter, be aware of the XY
## tolerance of the feature dataset containing your geometric network
#As gages being further than that tolerance will bring up "001191: No flag found from the flag junction feature class."

with arcpy.da.SearchCursor(gages_sel, ["site_no"]) as cursor:
    for row in cursor:
        #Print gage ID
        print row[0]
        expr = '"site_no" = ' + "'%s'" %row[0]
        #Select gage one by one based on gage ID
        arcpy.SelectLayerByAttribute_management("gage_sel_lyr", "NEW_SELECTION", expr)
        #Select all features in geometric network upstream of gage
        arcpy.TraceGeometricNetwork_management(in_geometric_network= network, out_network_layer = "up_trace_lyr", in_flags = "gage_sel_lyr", in_trace_task_type= "TRACE_UPSTREAM")
        up_tr = "up_trace_lyr/NHDFlowline_Network_nocoast2"
        #Iterate through each reach upstream of gage
        with arcpy.da.SearchCursor(up_tr, ["ComID"]) as cursor_up:
            for row_up in cursor_up:
                #Check whether ComID corresponds to ComID in the dam database
                if row_up[0] in dams_comid:
                    ##Makes sure that this ComID is not already in the dictionary for that gage (in the case of two dams being on the same reach)
                    if row[0] in gage_up_dic:
                        if row_up[0] not in gage_up_dic.get(row[0]):
                            #Write ComID to dictionary with the selected gage's ID as key and the reach with a dam's ComID as a value
                            gage_up_dic[row[0]].append(row_up[0])
                    else:
                        gage_up_dic[row[0]].append(row_up[0])
        #arcpy.CopyFeatures_management("up_trace_lyr/NHDFlowline_Network", "F:\Miscellaneous\Julian_Streamsites_gages\Gage_Dam_Pairing_uptrace.gdb\\up_%s" %(row[0]))
        #Delete temporary trace layer
        arcpy.Delete_management("up_trace_lyr")
arcpy.Delete_management("gage_sel_lyr")

gage_up_dic_copy = gage_up_dic
#Write out a csv file with a row for each pair of gage-reach with a dam
uptab = "F:\Miscellaneous\Julian_Streamsites_gages\gage_up_dic_25_20181204.csv"

with open(uptab, "wb") as csv_file:
    writer = csv.writer(csv_file)
    #Write headers
    writer.writerow(["SOURCE_FEA", "ComID_Dam"])
    for key in gage_up_dic.keys():
        for value in gage_up_dic[key]:
            writer.writerow([str(key),value])

#Identify all downstream dams

#Continuous subnetworks are those that are:
#1. in same basin
#2. same nearest downstream dam


########################################################################################################################
# EXTRA STUFF BY MATHIS
########################################################################################################################
# #Dissolve network by HYBAS_ID
# arcpy.env.workspace = dcigdb
# [f.name for f in arcpy.ListFields(netsub)]
# arcpy.Dissolve_management(netsub, 'netsubdiss', dissolve_field='PFAF_ID', statistics_fields=[['LENGTH_GEO', 'SUM']],
#                           multi_part='MULTI_PART', unsplit_lines='DISSOLVE_LINES')
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
