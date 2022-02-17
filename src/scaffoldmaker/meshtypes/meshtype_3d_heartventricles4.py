"""
Generates 3-D mesh of left and right ventricles below base plane.
Variant using collapsed/wedge elements at septum junction.
"""

from __future__ import division

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumElementIdentifier, getMaximumNodeIdentifier
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.heart_terms import get_heart_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_3d_heartarterialvalve1 import MeshType_3d_heartarterialvalve1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import createEllipsePoints, createEllipsoidPoints, getApproximateEllipsePerimeter, \
    getEllipseArcLength, getEllipseRadiansToX, createOvalPoints, revolvePoints
from scaffoldmaker.utils.interpolation import computeCubicHermiteDerivativeScaling, getCubicHermiteArcLength, \
    interpolateSampleCubicHermite, sampleCubicHermiteCurves, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.shieldmesh import ShieldMesh2D
from scaffoldmaker.utils.tracksurface import TrackSurface, calculate_surface_axes


class MeshType_3d_heartventricles4(Scaffold_base):
    '''
    Generates 3-D mesh of left and right ventricles below base plane.
    '''

    aorticValveScaffoldPackages = {
        'Default': ScaffoldPackage(MeshType_3d_heartarterialvalve1, {
            "rotation": [
                -42.280095308724924,
                20.287890155624396,
                -46.22275162169128
            ],
            "scaffoldSettings": {
                "Aortic": True,
                "Unit scale": 1.0
            },
            "scale": [
                0.32081136767341234,
                0.32081136767341234,
                0.32081136767341234
            ],
            "translation": [
                0.3053988636154118,
                0.09929523294480957,
                0.35979971157497675
            ]
        })
    }

    pulmonaryValveScaffoldPackages = {
        'Default': ScaffoldPackage(MeshType_3d_heartarterialvalve1, {
            "rotation": [
                -104.28829225942187,
                -14.895463238375616,
                39.25247517764923
            ],
            "scaffoldSettings": {
                "Aortic": False,
                "Unit scale": 1.0
            },
            "scale": [
                0.33536643732304455,
                0.33536643732304455,
                0.33536643732304455
            ],
            "translation": [
                0.3942377920137325,
                0.5283063303640532,
                0.24241414097833633
            ]
        })
    }

    @classmethod
    def getName(cls):
        return '3D Heart Ventricles 4'

    @classmethod
    def getParameterSetNames(cls):
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Unit Human 1',
            'Unit Mouse 1',
            'Unit Pig 1',
            'Unit Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        isHuman = 'Human' in parameterSetName
        isMouse = 'Mouse' in parameterSetName
        isPig = 'Pig' in parameterSetName
        isRat = 'Rat' in parameterSetName
        unitScale = 'Unit' in parameterSetName
        options = {}
        # first 3 are all even or all odd
        options['Number of elements across septum'] = 6  # S
        options['Number of elements around LV free wall'] = 8  # L
        options['Number of elements around RV free wall'] = 8  # R
        #options['Number of elements through LV wall'] = 1
        options['Number of elements up LV free wall'] = 5  # U
        options['Number of elements up LV apex'] = 0  # A each contributes 2 to number around LV free wall
        # N = L - 2A = 8
        # M = 2U - N = 6 (require 
        options['Unit scale'] = 1.0
        aorticValveOption = cls.aorticValveScaffoldPackages['Default']
        options['Aortic valve'] = copy.deepcopy(aorticValveOption)
        pulmonaryValveOption = cls.pulmonaryValveScaffoldPackages['Default']
        options['Pulmonary valve'] = copy.deepcopy(pulmonaryValveOption)
        options['Interventricular septum thickness'] = 0.12
        options['LV axis'] = [-0.2, 0.0, 1.0]
        options['LV axis apex height'] = 1.0
        options['LV axis base height'] = 0.5
        options['LV axis radius'] = 0.5
        options['LV apex thickness'] = 0.15
        options['LV free wall thickness'] = 0.15
        options['RV axis origin'] = [0.6, 0.0, 0.0]
        options['RV axis'] = [0.5, 0.1, 1.0]
        options['RV axis apex height'] = 1.0
        options['RV axis base height'] = 0.5
        options['RV axis radius'] = 0.5
        options['RV free wall thickness'] = 0.06
        options['RV shear xz'] = 0.0
        options['RV shear yz'] = 0.0
        options['RV taper'] = 0.0
        options['RV tilt degrees'] = 0.0

        #options['Use cross derivatives'] = False  # Removed from interface until working
        options['Refine'] = False
        options['Refine number of elements surface'] = 4
        options['Refine number of elements through wall'] = 1

        if isHuman:
            if not unitScale:
                options['Unit scale'] = 80.0
        elif isMouse or isRat:
            if not unitScale:
                options['Unit scale'] = 5.0 if isMouse else 12.0
        elif isPig:
            options['Number of elements up LV apex'] = 1
            if not unitScale:
                options['Unit scale'] = 80.0
            options['LV apex thickness'] = 0.07
            options['LV free wall thickness'] = 0.17
            options['Interventricular septum thickness'] = 0.14
            options['RV free wall thickness'] = 0.06
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            'Number of elements across septum',
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            #'Number of elements through LV wall',
            'Number of elements up LV free wall',
            'Number of elements up LV apex',
            'Unit scale',
            'Aortic valve',
            'Pulmonary valve',
            'Interventricular septum thickness',
            'LV axis',
            'LV axis apex height',
            'LV axis base height',
            'LV axis radius',
            'LV apex thickness',
            'LV free wall thickness',
            'RV axis origin',
            'RV axis',
            'RV axis apex height',
            'RV axis base height',
            'RV axis radius',
            'RV free wall thickness',
            'RV shear xz',
            'RV shear yz',
            'RV taper',
            'RV tilt degrees',
            # 'Use cross derivatives',  # Removed from interface until working
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName in ('Aortic valve', 'Pulmonary valve'):
            return [MeshType_3d_heartarterialvalve1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Aortic valve':
            return list(cls.aorticValveScaffoldPackages.keys())
        if optionName == 'Pulmonary valve':
            return list(cls.pulmonaryValveScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Aortic valve':
            if not parameterSetName:
                parameterSetName = list(cls.aorticValveScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.aorticValveScaffoldPackages[parameterSetName])
        if optionName == 'Pulmonary valve':
            if not parameterSetName:
                parameterSetName = list(cls.pulmonaryValveScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.pulmonaryValveScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def getElementsCounts(cls, options):
        elementsCountAcrossIVSeptum = options['Number of elements across septum']
        elementsCountAroundLVFreeWall = options['Number of elements around LV free wall']
        elementsCountAroundRVFreeWall = options['Number of elements around RV free wall']
        #elementsCountThroughLVWall = options['Number of elements through LV wall']
        elementsCountUpLVFreeWall = options['Number of elements up LV free wall']
        elementsCountUpLVApex = options['Number of elements up LV apex']
        nas = elementsCountAcrossIVSeptum//2
        nal = elementsCountAroundLVFreeWall//2 - elementsCountUpLVApex
        nar = elementsCountAroundRVFreeWall//2
        nul = elementsCountUpLVFreeWall
        # number around half = must be matched on RV and septum
        nah = nal + elementsCountUpLVFreeWall - 2
        elementsCountAroundFull = 2*nah + (options['Number of elements across septum'] % 2)
        # number up rv, min 2
        elementsCountUpRVFreeWall = nah + 2 - nar
        elementsCountUpIVSeptum = nah + 2 - elementsCountAcrossIVSeptum//2
        return elementsCountAroundFull, elementsCountUpLVApex, \
            elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
            elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
            elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall

    @classmethod
    def checkOptions(cls, options):
        '''
        :return:  True if dependent options changed, otherwise False.
        '''
        dependentChanges = False
        if options['Number of elements up LV apex'] < 0:
            options['Number of elements up LV apex'] = 0
        if options['Number of elements up LV free wall'] < (options['Number of elements up LV apex'] + 2):
            options['Number of elements up LV free wall'] = options['Number of elements up LV apex'] + 2
            dependentChanges = True
        for key in [
            'Number of elements across septum',
            'Number of elements around LV free wall',
            'Number of elements around RV free wall',
            ]:
            if options[key] < 4:
                options[key] = 4
        if options['Number of elements across septum'] % 2:
            if 0 == (options['Number of elements around LV free wall'] % 2):
                options['Number of elements around LV free wall'] += 1
                dependentChanges = True
            if 0 == (options['Number of elements around RV free wall'] % 2):
                options['Number of elements around RV free wall'] += 1
                dependentChanges = True
        nas = options['Number of elements across septum']//2
        nal = options['Number of elements around LV free wall']//2 - options['Number of elements up LV apex']
        nar = options['Number of elements around RV free wall']//2
        nul = options['Number of elements up LV free wall']
        # number around half = must be matched on RV and septum
        nah = nal + nul - 2
        # number up rv, min 2
        nur = nah + 2 - nar
        if nur < 2:
            options['Number of elements around RV free wall'] += 2*(2 - nur)
            dependentChanges = True
        # number up septum, min 2
        nus = nah + 2 - nas
        if nus < 2:
            options['Number of elements across septum'] += 2*(2 - nus)
            dependentChanges = True
        for key in [
            #'Number of elements through LV wall',
            'Refine number of elements surface',
            'Refine number of elements through wall',
            ]:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Unit scale',
            'Interventricular septum thickness',
            'LV apex thickness',
            'LV axis apex height',
            'LV axis base height',
            'LV axis radius',
            'LV free wall thickness',
            'RV axis apex height',
            'RV axis base height',
            'RV axis radius',
            'RV free wall thickness'
            ]:
            if options[key] < 0.0:
                options[key] = 0.0
        # for key in [
        #     ]:
        #     if options[key] < 0.01:
        #         options[key] = 0.01
        #     elif options[key] > 0.99:
        #         options[key] = 0.99
        for key in ['LV axis',
                    'RV axis',
                    'RV axis origin']:
            if len(options[key]) != 3:
                if len(options[key]) > 3:
                    options[key] = options[key][:3]
                else:
                    while len(options[key]) != 3:
                        options[key].append(0.0)
        for key in ['LV axis',
                    'RV axis']:
            mag = vector.magnitude(options[key])
            if mag < 0.00001:
                options[key] = [0.0, 0.0, 1.0]
        if options['RV tilt degrees'] < 0.0:
            options['RV tilt degrees'] = 0.0
        elif options['RV tilt degrees'] > 44.999:
            options['RV tilt degrees'] = 44.999
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        """
        elementsCountAroundFull, elementsCountUpLVApex, \
            elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
            elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
            elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall = cls.getElementsCounts(options)
        unitScale = options['Unit scale']
        aorticValve = options['Aortic valve']
        pulmonaryValve = options['Pulmonary valve']
        ivSeptumThickness = unitScale*options['Interventricular septum thickness']
        lvApexThickness = unitScale*options['LV apex thickness']
        lvAxis = vector.normalise(options['LV axis'])
        lvAxisApexHeight = unitScale*options['LV axis apex height']
        lvAxisBaseHeight = unitScale*options['LV axis base height']
        lvAxisRadius = unitScale*options['LV axis radius']
        lvFreeWallThickness = unitScale*options['LV free wall thickness']
        rvOrigin = [unitScale*s for s in options['RV axis origin']]
        rvAxis = vector.normalise(options['RV axis'])
        rvAxisApexHeight = unitScale*options['RV axis apex height']
        rvAxisBaseHeight = unitScale*options['RV axis base height']
        rvAxisRadius = unitScale*options['RV axis radius']
        rvFreeWallThickness = unitScale*options['RV free wall thickness']
        rvShearXZ = options['RV shear xz']
        rvShearYZ = options['RV shear yz']
        rvTaper = options['RV taper']
        rvTiltRadians = math.radians(options['RV tilt degrees'])

        useCrossDerivatives = False  # options['Use cross derivatives']  # Removed from interface until working

        #print("elementsCountAroundFull", elementsCountAroundFull)
        #print("elementsCountUpRVFreeWall", elementsCountUpRVFreeWall)

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        cache = fieldmodule.createFieldcache()

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldmodule.findMeshByDimension(3)

        heartGroup = AnnotationGroup(region, get_heart_term("heart"))
        apexGroup = AnnotationGroup(region, get_heart_term("apex of heart"))
        lvGroup = AnnotationGroup(region, get_heart_term("left ventricle myocardium"))
        rvGroup = AnnotationGroup(region, get_heart_term("right ventricle myocardium"))
        vSeptumGroup = AnnotationGroup(region, get_heart_term("interventricular septum"))
        annotationGroups = [ heartGroup, apexGroup, lvGroup, rvGroup, vSeptumGroup ]

        # annotation fiducial points
        markerGroup = findOrCreateFieldGroup(fieldmodule, "marker")
        markerName = findOrCreateFieldStoredString(fieldmodule, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fieldmodule, mesh, name="marker_location")

        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        #################
        # Create geometry
        #################

        aorticValve.generate(region, applyTransformationGroupName="root of aorta")
        annotationGroups += aorticValve.getAnnotationGroups()

        pulmonaryValve.generate(region, applyTransformationGroupName="root of pulmonary trunk")
        annotationGroups += pulmonaryValve.getAnnotationGroups()

        elementsCountAroundLVTrackSurface = 16
        elementsCountUpLVTrackSurface = 8

        # useHeight = min(max(0.0, lvOuterHeight), 2.0*lvAxisApexHeight)
        # baseRadiansUp = getEllipseRadiansToX(lvAxisApexHeight, 0.0, lvAxisApexHeight - useHeight,
        #                                      initialTheta=0.5*math.pi*useHeight/lvAxisApexHeight)
        # baseProportionUp = 2.0*getEllipseArcLength(lvAxisApexHeight, lvAxisRadius, 0.0, baseRadiansUp) \
        #                    / getApproximateEllipsePerimeter(lvAxisApexHeight, lvAxisRadius)
        lvOrigin = [0.0, 0.0, 0.0]
        lvApexAxis = vector.setMagnitude(lvAxis, -lvAxisApexHeight)
        lvSideAxis = vector.setMagnitude(vector.crossproduct3(lvAxis, [0.0, -1.0, 0.0]), lvAxisRadius)
        lx, ld = createOvalPoints(lvOrigin, lvApexAxis, lvSideAxis, lvAxisBaseHeight,
                                  elementsCountUpLVTrackSurface, 0.0, math.pi)

        ltx, ltd1, ltd2 = revolvePoints(lx, ld, lvOrigin, lvAxis, elementsCountAroundLVTrackSurface, loop=True)
        lvTrackSurface = TrackSurface(elementsCountAroundLVTrackSurface, elementsCountUpLVTrackSurface, ltx, ltd1, ltd2,
                                      loop1=True)

        elementsCountAroundRVTrackSurface = 16
        elementsCountUpRVTrackSurface = 8

        rvox, rvoy, startRadians, endRadians = \
            determineRVTiltOval(rvAxisApexHeight, rvAxisBaseHeight, rvAxisRadius, rvTiltRadians)

        rvPost = vector.normalise(vector.crossproduct3(rvAxis, [-1.0, 0.0, 0.0]))
        rvSide = vector.crossproduct3(rvPost, rvAxis)
        rvConstructionApex = [rvOrigin[c] - rvox*rvAxis[c] - rvoy*rvSide[c] for c in range(3)]
        sin_alpha = math.sin(rvTiltRadians)
        cos_alpha = math.cos(rvTiltRadians)
        rvApexAxis = vector.setMagnitude([cos_alpha*rvAxis[c] + sin_alpha*rvSide[c] for c in range(3)],
                                         -rvAxisApexHeight)
        scaledRvOuterRadius = (1.0 - rvTaper)*rvAxisRadius
        rvSideAxis = vector.setMagnitude(vector.crossproduct3(rvApexAxis, rvPost), scaledRvOuterRadius)
        rx, rd = createOvalPoints(rvOrigin, rvApexAxis, rvSideAxis, rvAxisBaseHeight,
                                  elementsCountUpRVTrackSurface, startRadians, endRadians)

        # taper only implemented for vertical RV axis
        scaledRvTaper = rvTaper/(1.0 - rvTaper)/rvAxisApexHeight
        for n in range(len(rx)):
            x = rx[n]
            dx = x[0] - rvConstructionApex[0]
            dz = x[2] - rvConstructionApex[2]
            x[0] += scaledRvTaper*dx*dz
            d = rd[n]
            d[0] += scaledRvTaper*(dx*d[2] + dz*d[0])

        rtx, rtd1, rtd2 = revolvePoints(rx, rd, rvConstructionApex, rvAxis, elementsCountAroundRVTrackSurface)

        # shear in xz, yz planes
        for n in range((elementsCountAroundRVTrackSurface + 1)*(elementsCountUpRVTrackSurface + 1)):
            x = rtx[n]
            dz = x[2] - rvConstructionApex[2]
            x[0] += rvShearXZ*dz
            x[1] += rvShearYZ*dz
            for d in (rtd1[n], rtd2[n]):
                d[0] += rvShearXZ*d[2]
                d[1] += rvShearYZ*d[2]

        rvTrackSurface = TrackSurface(elementsCountAroundRVTrackSurface, elementsCountUpRVTrackSurface, rtx, rtd1, rtd2)

        #################
        # Create nodes
        #################

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        drawLVTrackSurface = True
        drawRVTrackSurface = True

        if drawLVTrackSurface or drawRVTrackSurface:
            nodetemplate12 = nodes.createNodetemplate()
            nodetemplate12.defineField(coordinates)
            nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate12.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

        if drawLVTrackSurface:
            lvTrackSurfaceFirstNodeIdentifier = nodeIdentifier
            for n in range(len(ltx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate12)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, ltx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, ltd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, ltd2[n])
                nodeIdentifier += 1

        if drawRVTrackSurface:
            rvTrackSurfaceFirstNodeIdentifier = nodeIdentifier
            for n in range(len(rtx)):
                node = nodes.createNode(nodeIdentifier, nodetemplate12)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rtx [n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rtd1[n])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rtd2[n])
                nodeIdentifier += 1

        #################
        # Create elements
        #################

        heartMeshGroup = heartGroup.getMeshGroup(mesh)
        lvMeshGroup = lvGroup.getMeshGroup(mesh)
        rvMeshGroup = rvGroup.getMeshGroup(mesh)
        vSeptumMeshGroup = vSeptumGroup.getMeshGroup(mesh)

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        # elementIdentifier = lvShield.generateElements(fieldmodule, coordinates, elementIdentifier, [ heartMeshGroup, lvMeshGroup ])
        # elementIdentifier = rvShield.generateElements(fieldmodule, coordinates, elementIdentifier, [ heartMeshGroup, rvMeshGroup ])

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        if drawLVTrackSurface or drawRVTrackSurface:
            mesh2d = fieldmodule.findMeshByDimension(2)
            bicubicHermiteBasis = fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
            eft2d = mesh2d.createElementfieldtemplate(bicubicHermiteBasis)
            # remove cross derivative 12
            for n in range(4):
                r = eft2d.setFunctionNumberOfTerms(n*4 + 4, 0)
            elementtemplate2d = mesh2d.createElementtemplate()
            elementtemplate2d.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            elementtemplate2d.defineField(coordinates, -1, eft2d)

        if drawLVTrackSurface:
            nodesCount1 = elementsCountAroundLVTrackSurface
            for e2 in range(elementsCountUpLVTrackSurface):
                for e1 in range(elementsCountAroundLVTrackSurface):
                    element = mesh2d.createElement(-1, elementtemplate2d)  # since on 2-D mesh
                    nid1 = lvTrackSurfaceFirstNodeIdentifier + e2*nodesCount1 + e1
                    pid1 = lvTrackSurfaceFirstNodeIdentifier + e2*nodesCount1 + \
                           (e1 + 1) % elementsCountAroundLVTrackSurface
                    element.setNodesByIdentifier(eft2d, [nid1, pid1, nid1 + nodesCount1, pid1 + nodesCount1])

        if drawRVTrackSurface:
            nodesCount1 = elementsCountAroundRVTrackSurface + 1
            for e2 in range(elementsCountUpRVTrackSurface):
                for e1 in range(elementsCountAroundRVTrackSurface):
                    element = mesh2d.createElement(-1, elementtemplate2d)  # since on 2-D mesh
                    nid1 = rvTrackSurfaceFirstNodeIdentifier + e2*nodesCount1 + e1
                    element.setNodesByIdentifier(eft2d, [nid1, nid1 + 1, nid1 + nodesCount1, nid1 + nodesCount1 + 1])

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        Stops at end of ventricles, hence can be called from ventriclesbase.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        # assert isinstance(meshrefinement, MeshRefinement)
        # annotationGroups = meshrefinement.getAnnotationGroups()
        # elementsCountAroundFull, elementsCountUpLVApex, \
        #     elementsCountAcrossIVSeptum, elementsCountUpIVSeptum, \
        #     elementsCountAroundLVFreeWall, elementsCountUpLVFreeWall, \
        #     elementsCountAroundRVFreeWall, elementsCountUpRVFreeWall = cls.getElementsCounts(options)
        # refineElementsCountSurface = options['Refine number of elements surface']
        # #refineElementsCountThroughLVWall = options['Refine number of elements through LV wall']
        # refineElementsCountThroughWall = options['Refine number of elements through wall']
        # lvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left ventricle myocardium"))
        # rvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right ventricle myocardium"))
        # vSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interventricular septum"))
        # lvMeshGroup = lvGroup.getMeshGroup(meshrefinement._sourceMesh)
        # rvMeshGroup = rvGroup.getMeshGroup(meshrefinement._sourceMesh)
        # vSeptumMeshGroup = vSeptumGroup.getMeshGroup(meshrefinement._sourceMesh)
        # elementsCountVentricles = refineElementsCountThroughWall*(
        #     (elementsCountUpLVFreeWall - elementsCountUpLVApex - 2)*elementsCountAroundLVFreeWall + (elementsCountUpLVApex + 2)*(elementsCountAroundLVFreeWall - 2) +
        #     (elementsCountUpRVFreeWall - 2)*elementsCountAroundRVFreeWall + 2*(elementsCountAroundRVFreeWall - 2) +
        #     (elementsCountUpIVSeptum - 2)*elementsCountAcrossIVSeptum + 2*(elementsCountAcrossIVSeptum - 2) +
        #     elementsCountAroundFull)
        # element = meshrefinement._sourceElementiterator.next()
        # lastVentriclesElementIdentifier = element.getIdentifier() + elementsCountVentricles - 1
        # while element.isValid():
        #     elementIdentifier = element.getIdentifier()
        #     numberInXi1 = refineElementsCountSurface
        #     numberInXi2 = refineElementsCountSurface
        #     numberInXi3 = refineElementsCountThroughWall
        #     #if lvMeshGroup.containsElement(element):
        #     #    numberInXi3 = refineElementsCountThroughLVWall
        #     meshrefinement.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
        #     if elementIdentifier == lastVentriclesElementIdentifier:
        #         return  # finish on last so can continue in ventriclesbase
        #     element = meshrefinement._sourceElementiterator.next()
        pass


    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        # create endocardium and epicardium groups
        # fm = region.getFieldmodule()
        # lvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("left ventricle myocardium"))
        # rvGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("right ventricle myocardium"))
        # vSeptumGroup = getAnnotationGroupForTerm(annotationGroups, get_heart_term("interventricular septum"))
        # mesh2d = fm.findMeshByDimension(2)
        # is_exterior = fm.createFieldIsExterior()
        # is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        # is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        # is_lv = lvGroup.getFieldElementGroup(mesh2d)
        # is_rv = rvGroup.getFieldElementGroup(mesh2d)
        # is_lv_endo = fm.createFieldAnd(is_lv, is_exterior_face_xi3_0)
        # is_rv_endo = fm.createFieldOr(fm.createFieldAnd(fm.createFieldAnd(is_rv, is_exterior_face_xi3_0),
        #                                                 fm.createFieldNot(is_lv_endo)),
        #                               fm.createFieldAnd(vSeptumGroup.getFieldElementGroup(mesh2d), is_exterior_face_xi3_1))
        # is_v_epi = fm.createFieldAnd(fm.createFieldOr(is_lv, is_rv),
        #                              fm.createFieldAnd(is_exterior_face_xi3_1,
        #                                              fm.createFieldNot(vSeptumGroup.getFieldElementGroup(mesh2d))))
        # epiGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("epicardium"))
        # epiGroup.getMeshGroup(mesh2d).addElementsConditional(is_v_epi)
        # lvEndoGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("endocardium of left ventricle"))
        # lvEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_lv_endo)
        # rvEndoGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_heart_term("endocardium of right ventricle"))
        # rvEndoGroup.getMeshGroup(mesh2d).addElementsConditional(is_rv_endo)
        pass


def determineRVTiltOval(apex_length, base_length, side_radius, alpha):
    """
    Calculate oval origin relative to apex at (0,0) for tilt angle alpha,
    such that the apex point is tangential all around.
    :param apex_length: oval ellipse axis length down to apex
    :param base_length: ovel ellipse axis length up to base
    :param side_radius: ovel ellipse side radius
    :param alpha: Tilt angle in radians in [0, pi/4)
    :return: oval origin x, oval origin y, start angle, end angle.
    """
    assert 0.0 <= alpha < 0.25*math.pi
    if alpha == 0.0:
        return apex_length, 0.0, 0.0, math.pi
    # calculate theta = angle around ellipse to apex point
    a = apex_length
    b = side_radius
    sin_alpha = math.sin(alpha)
    cos_alpha = math.cos(alpha)
    tan_alpha = sin_alpha/cos_alpha
    theta = math.atan2(b*sin_alpha, a*cos_alpha)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    # tan_theta = sin_theta/cos_theta
    ay = b*sin_theta
    ax1 = a*(1.0 - cos_theta)
    ax2 = ay/tan_alpha
    aa = ay/sin_alpha
    ac = a - ax1 - ax2
    origin_x = ax2 + ac*cos_alpha
    origin_y = ac*sin_alpha
    start_angle = theta
    # GRC compute later:
    end_angle = 1.1*math.pi
    return origin_x, origin_y, start_angle, end_angle
