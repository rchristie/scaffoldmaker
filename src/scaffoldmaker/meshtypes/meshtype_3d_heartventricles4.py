"""
Generates 3-D mesh of left and right ventricles below base plane.
Variant using collapsed/wedge elements at septum junction.
"""

from __future__ import division

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
        options['Interventricular septum thickness'] = 0.12
        options['Interventricular sulcus apex transition length'] = 0.2
        options['Interventricular sulcus base transition length'] = 0.15
        options['LV apex thickness'] = 0.08
        options['LV axis apex height'] = 1.0
        options['LV axis base height'] = 0.5
        options['LV outer height'] = 1.0
        options['LV outer diameter'] = 1.0
        options['LV free wall thickness'] = 0.15
        options['RV apex relative coordinates'] = [0.6, -0.6, -0.8]
        options['RV free wall thickness'] = 0.05
        options['RV relative height'] = 0.8
        options['RV relative radius'] = 0.5
        options['RV rotation axis'] = [0.0, 0.0, 1.0]

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
            'Interventricular septum thickness',
            'Interventricular sulcus apex transition length',
            'Interventricular sulcus base transition length',
            'LV apex thickness',
            'LV axis apex height',
            'LV axis base height',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
            'RV apex relative coordinates',
            'RV free wall thickness',
            'RV relative height',
            'RV relative radius',
            'RV rotation axis',
            # 'Use cross derivatives',  # Removed from interface until working
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall'
        ]

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
            'Interventricular sulcus apex transition length',
            'Interventricular sulcus base transition length',
            'LV apex thickness',
            'LV axis apex height',
            'LV axis base height',
            'LV outer height',
            'LV outer diameter',
            'LV free wall thickness',
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
        for key in [
                'RV apex relative coordinates',
                'RV rotation axis']:
            if len(options[key]) != 3:
                if len(options[key]) > 3:
                    options[key] = options[key][:3]
                else:
                    while len(options[key]) != 3:
                        options[key].append(0.0)
        mag = vector.magnitude(options['RV rotation axis'])
        if mag < 0.00001:
            options['RV rotation axis'] = [0.0, 0.0, 1.0]
        elif not (0.99999 < mag < 1.00001):
            options['RV rotation axis'] = vector.setMagnitude(options['RV rotation axis'], 1.0)
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
        ivSeptumThickness = unitScale*options['Interventricular septum thickness']
        ivSulcusApexTransitionLength = options['Interventricular sulcus apex transition length']
        ivSulcusBaseTransitionLength = options['Interventricular sulcus base transition length']
        lvApexThickness = unitScale*options['LV apex thickness']
        lvAxisApexHeight = unitScale*options['LV axis apex height']
        lvAxisBaseHeight = unitScale*options['LV axis base height']
        lvOuterHeight = unitScale*options['LV outer height']
        lvOuterRadius = unitScale*0.5*options['LV outer diameter']
        lvFreeWallThickness = unitScale*options['LV free wall thickness']
        rvApexRelativeCoordinates = options['RV apex relative coordinates']
        rvFreeWallThickness = unitScale*options['RV free wall thickness']
        rvAxisHeight = lvAxisApexHeight * options['RV relative height']
        rvAxisWidth = lvOuterRadius * options['RV relative radius']
        rvRotationAxis = options['RV rotation axis']
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

        elementsCountAroundLVTrackSurface = 16
        elementsCountUpLVTrackSurface = 8

        # useHeight = min(max(0.0, lvOuterHeight), 2.0*lvAxisApexHeight)
        # baseRadiansUp = getEllipseRadiansToX(lvAxisApexHeight, 0.0, lvAxisApexHeight - useHeight,
        #                                      initialTheta=0.5*math.pi*useHeight/lvAxisApexHeight)
        # baseProportionUp = 2.0*getEllipseArcLength(lvAxisApexHeight, lvOuterRadius, 0.0, baseRadiansUp) \
        #                    / getApproximateEllipsePerimeter(lvAxisApexHeight, lvOuterRadius)
        origin = [0.0, 0.0, 0.0]
        lx, ld = createOvalPoints(origin, [0.0, 0.0, -lvAxisApexHeight], [0.0, -lvOuterRadius, 0.0], lvAxisBaseHeight,
                                  elementsCountUpLVTrackSurface, math.pi)

        ltx, ltd1, ltd2 = revolvePoints(lx, ld, origin, [0.0, 0.0, 1.0], elementsCountAroundLVTrackSurface,
                                        loop=True)

        lvTrackSurface = TrackSurface(elementsCountAroundLVTrackSurface, elementsCountUpLVTrackSurface, ltx, ltd1, ltd2,
                                      loop1=True)

        elementsCountAroundRVTrackSurface = 16
        elementsCountUpRVTrackSurface = 8

        rvApexCoordinates = [rvApexRelativeCoordinates[0] * lvOuterRadius,
                             rvApexRelativeCoordinates[1] * lvOuterRadius,
                             rvApexRelativeCoordinates[2] * lvAxisApexHeight]
        centre = [rvApexCoordinates[0], rvApexCoordinates[1], rvApexCoordinates[2] + rvAxisHeight]

        rx, rd = createEllipsePoints(centre, [0.0, 0.0, -rvAxisHeight], [0.0, -rvAxisWidth, 0.0],
                                     elementsCountUpRVTrackSurface, 1.25*math.pi)
        rtx, rtd1, rtd2 = revolvePoints(rx, rd, rvApexCoordinates, rvRotationAxis, elementsCountAroundRVTrackSurface, \
                                        1.25*math.pi, startRadians=-0.25*math.pi)

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
