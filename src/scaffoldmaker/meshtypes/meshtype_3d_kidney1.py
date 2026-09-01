"""
Generates a 3D kidney capsule.
"""
from cmlibs.maths.vectorops import add, cross, div, magnitude, mult, normalize, set_magnitude, sub
from cmlibs.utils.zinc.field import find_or_create_field_coordinates, findOrCreateFieldCoordinates
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm, findAnnotationGroupByName
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.ellipsoidmesh import EllipsoidMesh
from scaffoldmaker.utils.geometry import sampleEllipsePoints
from scaffoldmaker.utils.hextetrahedronmesh import HexTetrahedronMesh
from scaffoldmaker.utils.interpolation import (
    interpolateSampleCubicHermite, sampleHermiteCurve, smoothCubicHermiteDerivativesLine, sampleCubicHermiteCurves,
    smoothCurveSideCrossDerivatives)
from scaffoldmaker.utils.meshgeneratedata import MeshGenerateData
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.quadtrianglemesh import QuadTriangleMesh
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.zinc_utils import translate_nodeset_coordinates
import math


class KidneyTubeNetworkMeshGenerateData(TubeNetworkMeshGenerateData):

    def __init__(self, region, meshDimension, coordinateFieldName="coordinates",
                 startNodeIdentifier=1, startElementIdentifier=1, isLinearThroughShell=False, isShowTrimSurfaces=False):
        """
        :param isLinearThroughWall: Callers should only set if 3-D with no core.
        :param isShowTrimSurfaces: Tells junction generateMesh to make 2-D trim surfaces.
        """
        super(KidneyTubeNetworkMeshGenerateData, self).__init__(
            region, meshDimension, coordinateFieldName, startNodeIdentifier, startElementIdentifier,
            isLinearThroughShell, isShowTrimSurfaces)
        # force these names for standard annotation groups in the base class
        self._coreGroup = self.getOrCreateAnnotationGroup(get_kidney_term("renal medulla"))
        self._shellGroup = self.getOrCreateAnnotationGroup(get_kidney_term("cortex of kidney"))

    def getMedullaMeshGroup(self):
        return self._coreGroup.getMeshGroup(self._mesh)

    def getCortexMeshGroup(self):
        return self._shellGroup.getMeshGroup(self._mesh)


class KidneyTubeNetworkMeshBuilder(TubeNetworkMeshBuilder):
    """
    Specialization of TubeNetworkMeshBuilder adding annotations for anterior, posterior, lateral, medial, and hilum regions.
    """

    def __init__(self, networkMesh: NetworkMesh, targetElementDensityAlongLongestSegment: float,
                 layoutAnnotationGroups: list=[], annotationElementsCountsAlong: list=[],
                 defaultElementsCountAround: int=8, annotationElementsCountsAround: list=[],
                 elementsCountThroughShell: int=1, isCore=False, elementsCountTransition: int=1,
                 defaultElementsCountCoreBoxMinor: int=2, annotationElementsCountsCoreBoxMinor: list=[],
                 defaultCoreBoundaryScalingMode=1, annotationCoreBoundaryScalingMode=[],
                 useOuterTrimSurfaces=True, showKidneys=[]):
        """
        Builds specialized continuous tube network meshes for kidney scaffold.
        :param showKidneys: List of flags for showing left and/or right kidneys.
        """
        super(KidneyTubeNetworkMeshBuilder, self).__init__(
            networkMesh, targetElementDensityAlongLongestSegment, layoutAnnotationGroups, annotationElementsCountsAlong,
            defaultElementsCountAround, annotationElementsCountsAround, elementsCountThroughShell, isCore,
            elementsCountTransition, defaultElementsCountCoreBoxMinor, annotationElementsCountsCoreBoxMinor,
            defaultCoreBoundaryScalingMode, annotationCoreBoundaryScalingMode, useOuterTrimSurfaces)

        self._showKidneys = showKidneys


    def generateMesh(self, generateData):
        super(KidneyTubeNetworkMeshBuilder, self).generateMesh(generateData)
        # build anterior, posterior, lateral, medial annotation groups
        anteriorMeshGroup = generateData.getAnteriorMeshGroup()
        posteriorMeshGroup = generateData.getPosteriorMeshGroup()
        lateralMeshGroup = generateData.getLateralMeshGroup()
        medialMeshGroup = generateData.getMedialMeshGroup()
        dorsalMeshGroup = generateData.getDorsalMeshGroup()
        ventralMeshGroup = generateData.getVentralMeshGroup()
        openingMeshGroup = generateData.getOpeningMeshGroup()

        elementsCountAround = self._defaultElementsCountAround
        halfElementsCountAround = elementsCountAround // 2
        increment = max(1, elementsCountAround // 8)

        leftKidney, rightKidney = 0, 1
        # Kidney configuration mapping
        kidney_configs = {
            leftKidney: {
                'lateral_flag': False,
                'medial_flag': True,
                'e1_start': halfElementsCountAround - increment,
                'e1_end': halfElementsCountAround + increment
            },
            rightKidney: {
                'lateral_flag': True,
                'medial_flag': False,
                'e1_start': -increment,
                'e1_end': increment
            }
        }

        def add_common_elements(mesh_obj, method_prefix=""):
            """
            Add D1 and D3 elements that are common to all kidneys.
            :param mesh_obj: The mesh object (segment or capMesh) to add elements to
            :param method_prefix: Prefix for method names (e.g., "addCap" for capMesh methods)
            """
            d1_method = f"{method_prefix}SideD1ElementsToMeshGroup"
            d3_method = f"{method_prefix}SideD3ElementsToMeshGroup"

            getattr(mesh_obj, d1_method)(False, anteriorMeshGroup)
            getattr(mesh_obj, d1_method)(True, posteriorMeshGroup)
            getattr(mesh_obj, d3_method)(False, ventralMeshGroup)
            getattr(mesh_obj, d3_method)(True, dorsalMeshGroup)

        def add_kidney_specific_elements(mesh_obj, kidney, method_prefix=""):
            """
            Add kidney-specific D2 elements.
            :param mesh_obj: The mesh object (segment or capMesh) to add elements to
            :param kidney: Kidney identifier (leftKidney=0, rightKidney=1)
            :param method_prefix: Prefix for method names (e.g., "addCap" for capMesh methods)
            """
            if kidney not in kidney_configs:
                return

            config = kidney_configs[kidney]
            d2_method = f"{method_prefix}SideD2ElementsToMeshGroup"

            getattr(mesh_obj, d2_method)(config['lateral_flag'], lateralMeshGroup)
            getattr(mesh_obj, d2_method)(config['medial_flag'], medialMeshGroup)

            # Shell opening elements only for segment (not capMesh)
            if "Cap" not in method_prefix:
                mesh_obj.addShellOpeningElementsToMeshGroup(config['e1_start'], config['e1_end'], openingMeshGroup)

        for kidney in [leftKidney, rightKidney]:
            if not self._showKidneys[kidney]:
                continue

            idx = 0 if False in self._showKidneys else kidney
            networkSegment = self._networkMesh.getNetworkSegments()[idx]
            segment = self._segments[networkSegment]
            segmentCaps = segment.getIsCap()
            capMesh = segment.getCapMesh() if True in segmentCaps else None

            # Apply common elements to segment
            add_common_elements(segment, "add")

            # Apply kidney-specific elements to segment
            add_kidney_specific_elements(segment, kidney, "add")

            # Apply elements to capMesh if it exists
            if capMesh:
                add_common_elements(capMesh, "addCap")
                add_kidney_specific_elements(capMesh, kidney, "addCap")



# class MeshType_1d_kidney_network_layout1(MeshType_1d_network_layout1):
#     """
#     Defines kidney network layout.
#     """
#
#     showKidneys = [False, False]
#
#     @classmethod
#     def getName(cls):
#         return "1D Kidney Network Layout 1"
#
#     @classmethod
#     def getParameterSetNames(cls):
#         return ["Default"]
#
#     @classmethod
#     def getDefaultOptions(cls, parameterSetName="Default"):
#         options = {}
#         options["Base parameter set"] = "Human 1" if (parameterSetName == "Default") else parameterSetName
#         options["Define inner coordinates"] = True
#         options["Left kidney"] = True
#         options["Right kidney"] = True
#         options["Kidney length"] = 1.0
#         options["Kidney width"] = 0.5
#         options["Kidney thickness"] = 0.4
#         options["Medulla to cortex proportion"] = 0.6
#         return options
#
#     @classmethod
#     def getOrderedOptionNames(cls):
#         return [
#             "Left kidney",
#             "Right kidney",
#             "Kidney length",
#             "Kidney width",
#             "Kidney thickness",
#             "Medulla to cortex proportion"
#         ]
#
#     @classmethod
#     def checkOptions(cls, options):
#         dependentChanges = False
#         for key in [
#             "Kidney length",
#             "Kidney width",
#             "Kidney thickness"
#         ]:
#             if options[key] < 0.1:
#                 options[key] = 0.1
#
#         if options["Medulla to cortex proportion"] < 0.1:
#             options["Medulla to cortex proportion"] = 0.1
#         elif options["Medulla to cortex proportion"] > 0.9:
#             options["Medulla to cortex proportion"] = 0.9
#
#         if not options["Left kidney"] and not options["Right kidney"]:
#             dependentChanges = True
#             options["Left kidney"] = True
#
#         return dependentChanges
#
#     @classmethod
#     def generateBaseMesh(cls, region, options):
#         """
#         Generate the unrefined mesh.
#         :param region: Zinc region to define model in. Must be empty.
#         :param options: Dict containing options. See getDefaultOptions().
#         :return: [] empty list of AnnotationGroup, NetworkMesh
#         """
#         # parameters
#         structure = options["Structure"] = cls.getLayoutStructure(options)
#         isLeftKidney = options["Left kidney"]
#         isRightKidney = options["Right kidney"]
#         kidneyLength = options["Kidney length"]
#         halfKidneyLength = 0.5 * kidneyLength
#         halfKidneyWidth = 0.5 * options["Kidney width"]
#         halfKidneyThickness = 0.5 * options["Kidney thickness"]
#         innerProportionDefault = options["Medulla to cortex proportion"]
#         cls.setShowKidneys(options)
#
#         networkMesh = NetworkMesh(structure)
#         networkMesh.create1DLayoutMesh(region)
#
#         fieldmodule = region.getFieldmodule()
#         mesh = fieldmodule.findMeshByDimension(1)
#
#         # set up element annotations
#         kidneyGroup = AnnotationGroup(region, get_kidney_term("kidney"))
#         kidneyMeshGroup = kidneyGroup.getMeshGroup(mesh)
#
#         leftKidneyGroup = AnnotationGroup(region, get_kidney_term("left kidney"))
#         leftKidneyMeshGroup = leftKidneyGroup.getMeshGroup(mesh)
#
#         rightKidneyGroup = AnnotationGroup(region, get_kidney_term("right kidney"))
#         rightKidneyMeshGroup = rightKidneyGroup.getMeshGroup(mesh)
#
#         annotationGroups = [kidneyGroup, leftKidneyGroup, rightKidneyGroup]
#         meshGroups = [kidneyMeshGroup, leftKidneyMeshGroup, rightKidneyMeshGroup]
#
#         # set coordinates (outer)
#         fieldcache = fieldmodule.createFieldcache()
#         coordinates = find_or_create_field_coordinates(fieldmodule)
#         # need to ensure inner coordinates are at least defined:
#         cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=innerProportionDefault)
#         innerCoordinates = find_or_create_field_coordinates(fieldmodule, "inner coordinates")
#         nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
#
#         # Kidney
#         nodeIdentifier = 1
#         elementIdentifier = 1
#         kidneyElementsCount = 2
#         capRadius = cls.getCapRadius(halfKidneyWidth, halfKidneyThickness) * (
#                 halfKidneyWidth * 0.45 + halfKidneyThickness * 0.55)
#         extensionLength = 0.5 * (halfKidneyWidth * 0.45 + halfKidneyThickness * 0.55)
#         halfLayoutLength = (halfKidneyLength - capRadius - extensionLength)
#         kidneyScale = 2 * halfLayoutLength / kidneyElementsCount
#
#         leftKidney, rightKidney = 0, 1
#         kidneys = [kidney for show, kidney in [(isLeftKidney, leftKidney), (isRightKidney, rightKidney)] if show]
#         for kidney in kidneys:
#             mx = [0.0, 0.0, 0.0]
#             d1 = [kidneyScale, 0.0, 0.0]
#             d3 = [0.0, 0.0, halfKidneyThickness]
#             id3 = mult(d3, innerProportionDefault)
#
#             tx = halfLayoutLength
#             sx = [-tx, 0.0, 0.0] if kidney is leftKidney else [-tx, 0.0, 0.0]
#             ex = [tx, 0.0, 0.0] if kidney is leftKidney else [tx, 0.0, 0.0]
#             sd1 = mult([-1.0, 0.0, 0.0], kidneyScale)
#             ed1 = [-sd1[0], sd1[1], sd1[2]]
#             nx, nd1 = sampleCubicHermiteCurves([sx, mx, ex], [sd1, d1, ed1], kidneyElementsCount)[0:2]
#             nd1 = smoothCubicHermiteDerivativesLine(nx, nd1)
#
#             sd2_list = []
#             sd3_list = []
#             sNodeIdentifiers = []
#             for e in range(kidneyElementsCount + 1):
#                 sNodeIdentifiers.append(nodeIdentifier)
#                 node = nodes.findNodeByIdentifier(nodeIdentifier)
#                 fieldcache.setNode(node)
#                 sd2 = set_magnitude(cross(d3, nd1[e]), halfKidneyWidth)
#                 sid2 = mult(sd2, innerProportionDefault)
#                 sd2_list.append(sd2)
#                 sd3_list.append(d3)
#                 for field, derivatives in ((coordinates, (nd1[e], sd2, d3)), (innerCoordinates, (nd1[e], sid2, id3))):
#                     setNodeFieldParameters(field, fieldcache, nx[e], *derivatives)
#                 nodeIdentifier += 1
#
#             sd12 = smoothCurveSideCrossDerivatives(nx, nd1, [sd2_list])[0]
#             sd13 = smoothCurveSideCrossDerivatives(nx, nd1, [sd3_list])[0]
#             for e in range(kidneyElementsCount + 1):
#                 node = nodes.findNodeByIdentifier(sNodeIdentifiers[e])
#                 fieldcache.setNode(node)
#                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[e])
#                 coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[e])
#                 sid12 = mult(sd12[e], innerProportionDefault)
#                 sid13 = mult(sd13[e], innerProportionDefault)
#                 innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sid12)
#                 innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sid13)
#
#             # add annotations
#             for e in range(kidneyElementsCount):
#                 element = mesh.findElementByIdentifier(elementIdentifier)
#                 meshGroups[0].addElement(element)
#                 if kidney is leftKidney and isLeftKidney:
#                     meshGroups[1].addElement(element)
#                 if kidney is rightKidney and isRightKidney:
#                     meshGroups[2].addElement(element)
#                 elementIdentifier += 1
#
#         return annotationGroups, networkMesh
#
#     @classmethod
#     def getLayoutStructure(cls, options):
#         """
#         Generate 1D layout structure based on the number of elements count along.
#         :param options: Dict containing options. See getDefaultOptions().
#         :return: string version of the 1D layout structure
#         """
#         nodes_count = 3
#         assert nodes_count > 1
#
#         left = f"({'-'.join(map(str, range(1, nodes_count + 1)))})"
#         if options["Left kidney"] and options["Right kidney"]:
#             right = f"({'-'.join(map(str, range(nodes_count + 1, 2 * nodes_count + 1)))})"
#             return f"{left},{right}"
#         return left
#
#     @classmethod
#     def getCapRadius(cls, majorRadius, minorRadius):
#         """
#         Calculate the radius of the cap mesh based on the major radius and the minor radius of tube cross-section.
#         :param majorRadius: The radius of a tube in the major-axis.
#         :param minorRadius: The radius of a tube in the minor-axis.
#         :return: Cap radius
#         """
#         if majorRadius > minorRadius:
#             return math.pow((majorRadius / minorRadius), 1 / 3)
#         elif majorRadius < minorRadius:
#             return math.pow((minorRadius / majorRadius), 1 / 3)
#         else:
#             return majorRadius
#
#     @classmethod
#     def getShowKidneys(cls):
#         return cls.showKidneys
#
#     @classmethod
#     def setShowKidneys(cls, options):
#         cls.showKidneys[0] = True if options["Left kidney"] else False
#         cls.showKidneys[1] = True if options["Right kidney"] else False


class MeshType_3d_kidney1(Scaffold_base):
    """
    Generates a 3-D Kidney.
    """

    @classmethod
    def getName(cls):
        return "3D Kidney 1"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Base parameter set"] = useParameterSetName
        options["Left kidney"] = True
        options["Right kidney"] = True
        options["Number of elements around"] = 16
        options["Number of elements along"] = 12
        options["Number of elements across core box minor"] = 4
        options["Number of elements across core transition"] = 1
        options["Number of elements through cortex"] = 2
        options["Kidney curvature"] = 1.0
        options["Kidney length"] = 1.0
        options["Kidney spacing"] = 1.0
        options["Kidney thickness"] = 0.4
        options["Kidney width"] = 0.5
        options["Medulla to cortex proportion"] = 0.6

        options["Refine"] = False
        options["Refine number of elements"] = 4
        options["Refine number of elements through cortex"] = 2
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Left kidney",
            "Right kidney",
            "Number of elements around",
            "Number of elements along",
            "Number of elements across core box minor",
            "Number of elements across core transition",
            "Number of elements through cortex",
            "Kidney curvature",
            "Kidney length",
            "Kidney spacing",
            "Kidney thickness",
            "Kidney width",
            "Refine",
            "Refine number of elements",
            "Refine number of elements through cortex"
        ]
        return optionNames

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False

        if options["Number of elements around"] < 8:
            options["Number of elements around"] = 8
        elif options["Number of elements around"] % 4:
            options["Number of elements around"] += 4 - (options["Number of elements around"] % 4)

        if options["Number of elements along"] < 4:
            options["Number of elements along"] = 4

        maxElementsCountCoreBoxMinor = options["Number of elements around"] // 2 - 2
        if options["Number of elements across core box minor"] < 2:
            options["Number of elements across core box minor"] = 2
        elif options["Number of elements across core box minor"] > maxElementsCountCoreBoxMinor:
            options["Number of elements across core box minor"] = maxElementsCountCoreBoxMinor
            dependentChanges = True
        elif options["Number of elements across core box minor"] % 2:
            options["Number of elements across core box minor"] += 1

        for key in [
            "Number of elements across core transition",
            "Number of elements through cortex",
            "Refine number of elements",
            "Refine number of elements through cortex"
        ]:
            if options[key] < 1:
                options[key] = 1

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        # networkLayout = options["Kidney network layout"]
        # layoutRegion = region.createRegion()
        # networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        # layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        # networkMesh = networkLayout.getConstructionObject()
        # showKidneys = getShowKidneysSettings()
        # isLeftKidney = showKidneys[0]
        # isRightKidney = showKidneys[1]

        outer_a = 0.25
        outer_b = 0.2
        outer_c = 0.5
        cortex_thickness = 0.3 * outer_a
        # sinus_proportion = 0.45
        sinus_a = outer_a * 0.45
        sinus_b = outer_b * 0.4
        sinus_c = outer_c * 0.5
        hilum_angle = math.radians(135)

        element_counts = [8, 8, 12]
        shell_count = 0
        transition_count = 2
        nway_d_factor = 0.6
        ellipsoid = EllipsoidMesh(element_counts, shell_count, transition_count, core=True)

        octant = build_kidney_capsule_octant(
            outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle,
            element_counts, transition_count, nway_d_factor=nway_d_factor)

        ellipsoid.merge_octant_plus1_quadrant(octant, 0)

        node_identifier = 1
        element_identifier = 1

        generate_data = MeshGenerateData(region, 3, "coordinates", node_identifier, element_identifier,)
        ellipsoid.generate_nodes(generate_data)
        ellipsoid.generate_elements(generate_data)

        annotation_groups = []

        # kidneyTubeNetworkMeshBuilder = KidneyTubeNetworkMeshBuilder(
        #     networkMesh,
        #     targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
        #     defaultElementsCountAround=options["Number of elements around"],
        #     elementsCountThroughShell=options["Number of elements through cortex"],
        #     layoutAnnotationGroups=layoutAnnotationGroups,
        #     isCore=True,
        #     elementsCountTransition=options["Number of elements across core transition"],
        #     defaultElementsCountCoreBoxMinor=options["Number of elements across core box minor"],
        #     annotationElementsCountsCoreBoxMinor=options["Annotation numbers of elements across core box minor"],
        #     showKidneys=showKidneys
        # )
        #
        # kidneyTubeNetworkMeshBuilder.build()
        # generateData = TubeNetworkMeshGenerateData(
        #     region, 3,
        #     isLinearThroughShell=False)
        # kidneyTubeNetworkMeshBuilder.generateMesh(generateData)
        # annotationGroups = generateData.getAnnotationGroups()
        #
        # # add kidney-specific annotation groups
        # fm = region.getFieldmodule()
        # coordinates = findOrCreateFieldCoordinates(fm)
        # mesh = generateData.getMesh()
        # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        #
        # coreGroup = getAnnotationGroupForTerm(annotationGroups, ("core", "")).getGroup()
        # shellGroup = getAnnotationGroupForTerm(annotationGroups, ("shell", "")).getGroup()
        # openingGroup = getAnnotationGroupForTerm(annotationGroups, ("opening", "")).getGroup()
        #
        # kidneyGroup = AnnotationGroup(region, get_kidney_term("kidney"))
        # kidneyNodesetGroup = kidneyGroup.getNodesetGroup(nodes)
        #
        # leftKidneyGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("left kidney"))
        # leftKidneyNodesetGroup = leftKidneyGroup.getNodesetGroup(nodes)
        #
        # rightKidneyGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("right kidney"))
        # rightKidneyNodesetGroup = rightKidneyGroup.getNodesetGroup(nodes)
        #
        # hilumGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("hilum of kidney"))
        # hilumGroup.getMeshGroup(mesh).addElementsConditional(openingGroup)
        #
        # cortexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("cortex of kidney"))
        # tempGroup = fm.createFieldSubtract(shellGroup, openingGroup)
        # cortexGroup.getMeshGroup(mesh).addElementsConditional(tempGroup)
        #
        # medullaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("renal medulla"))
        # medullaGroup.getMeshGroup(mesh).addElementsConditional(coreGroup)
        #
        # for term in ["core", "shell", "opening"]:
        #     annotationGroups.remove(findAnnotationGroupByName(annotationGroups, term))
        #
        # # marker points
        # leftSuperiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("superior pole of left kidney"))
        # leftInferiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("inferior pole of left kidney"))
        #
        # rightSuperiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("superior pole of right kidney"))
        # rightInferiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("inferior pole of right kidney"))
        #
        # markerList = []
        # elementsCountAround = options["Number of elements around"]
        # elementsCountAlong = int(options["Target element density along longest segment"] + 2) # extra sections from the cap mesh
        # elementsCountThroughShell = options["Number of elements through cortex"]
        # elementsCountCoreBoxMinor = options["Number of elements across core box minor"]
        # elementsCountCoreBoxMajor = (elementsCountAround // 2) - elementsCountCoreBoxMinor
        # elementsCountTransition = options["Number of elements across core transition"]
        #
        # box_count = elementsCountCoreBoxMinor * elementsCountCoreBoxMajor
        # depth = elementsCountThroughShell + elementsCountTransition
        # offset = elementsCountCoreBoxMajor // 2 * elementsCountCoreBoxMinor + elementsCountCoreBoxMinor // 2 + 1
        # cap_count = box_count * (depth + 1) + elementsCountAround * depth
        # tube_section_count = box_count + elementsCountAround * depth
        # tube_count = tube_section_count * elementsCountAlong
        # kidney_elements_count = cap_count * 2 + tube_count if showKidneys[0] else 0
        #
        # if isLeftKidney:
        #     idx = box_count * depth + offset
        #     markerList.append({"group": leftSuperiorPoleGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})
        #
        #     idx = cap_count + tube_count + (box_count * depth + offset)
        #     markerList.append({"group": leftInferiorPoleGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})
        #
        # if isRightKidney:
        #     idx = kidney_elements_count + box_count * depth + offset
        #     markerList.append({"group": rightSuperiorPoleGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})
        #
        #     idx = kidney_elements_count + cap_count + tube_count + (box_count * depth + offset)
        #     markerList.append({"group": rightInferiorPoleGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})
        #
        # nodeIdentifier = generateData.nextNodeIdentifier()
        # for marker in markerList:
        #     annotationGroup = marker["group"]
        #     markerNode = annotationGroup.createMarkerNode(
        #         nodeIdentifier, element=mesh.findElementByIdentifier(marker["elementId"]), xi=marker["xi"])
        #     annotationGroup.setMarkerMaterialCoordinates(coordinates)
        #     kidneyNodesetGroup.addNode(markerNode)
        #     nodeIdentifier += 1
        #
        # # transformation
        # leftKidney, rightKidney = 0, 1
        # kidneys = [kidney for show, kidney in [(isLeftKidney, leftKidney), (isRightKidney, rightKidney)] if show]
        # for kidney in kidneys:
        #     isLeft = True if kidney == leftKidney else False
        #     isRight = True if kidney == rightKidney else False
        #     spacing = options["Kidney spacing"] / 2 if isLeft else -options["Kidney spacing"] / 2
        #     curvature = -options["Kidney curvature"] if isLeft else options["Kidney curvature"]
        #
        #     kidneyNodeset = leftKidneyNodesetGroup if isLeft else rightKidneyNodesetGroup
        #
        #     if curvature != 0.0:
        #         if isLeft:
        #             bendKidneyMeshAroundZAxis(curvature, fm, coordinates, kidneyNodeset, stationaryPointXY=[0.05, 0.0])
        #         if isRight:
        #             bendKidneyMeshAroundZAxis(curvature, fm, coordinates, kidneyNodeset, stationaryPointXY=[-0.05, 0.0])
        #
        #     translate_nodeset_coordinates(kidneyNodeset, coordinates, [0.0, spacing, 0.0])

        return annotation_groups, None


    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refine_count = options["Refine number of elements"]
        meshRefinement.refineAllElementsCubeStandard3d(refine_count, refine_count, refine_count)


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
        return
        show_kidneys = getShowKidneysSettings()

        # Initialize field module and meshes
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)
        mesh3d = fm.findMeshByDimension(3)
        is_exterior = fm.createFieldIsExterior()

        # Get base kidney group
        kidney_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("kidney")).getGroup()

        # Create side groups and tracking dictionaries
        side_groups = {}
        side_kidney_groups = {"left": {}, "right": {}}

        for side in ["left", "right"]:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(f"{side} kidney"))
            side_groups[side] = group.getGroup()
            side_kidney_groups[side][f"{side} kidney"] = group

        # Create kidney part groups (cortex, hilum, medulla)
        create_kidney_part_groups(fm, mesh3d, annotationGroups, region, side_groups, side_kidney_groups)

        # Create capsule groups
        create_capsule_groups(fm, mesh2d, annotationGroups, region, kidney_group, side_groups, side_kidney_groups, is_exterior)

        # Create surface and edge annotation groups
        create_surface_and_edge_groups(fm, mesh2d, mesh1d, annotationGroups, region, side_groups, side_kidney_groups, is_exterior)

        # Remove groups based on kidney visibility settings
        remove_hidden_kidney_groups(show_kidneys, side_kidney_groups, annotationGroups)


# def getShowKidneysSettings():
#     return MeshType_1d_kidney_network_layout1.getShowKidneys()


def create_kidney_part_groups(fm, mesh3d, annotationGroups, region, side_groups, side_kidney_groups):
    """
    Create cortex, hilum, and medulla groups for each kidney side.
    """
    kidney_parts = ["cortex", "hilum", "medulla"]

    for part in kidney_parts:
        # Get the anatomical term for the part
        arb_term = f"renal {part}" if part == "medulla" else f"{part} of kidney"
        arb_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term(arb_term)).getGroup()

        # Create side-specific part groups
        for side in ["left", "right"]:
            part_term = f"{part} of {side} kidney"
            part_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(part_term))
            part_group.getMeshGroup(mesh3d).addElementsConditional(fm.createFieldAnd(arb_group, side_groups[side]))
            side_kidney_groups[side][part_term] = part_group


def create_capsule_groups(fm, mesh2d, annotationGroups, region, kidney_group, side_groups, side_kidney_groups, is_exterior):
    """
    Create kidney capsule surface groups.
    """
    # General kidney capsule group
    kidney_capsule_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("kidney capsule"))
    kidney_exterior = fm.createFieldAnd(kidney_group, is_exterior)
    kidney_capsule_group.getMeshGroup(mesh2d).addElementsConditional(kidney_exterior)

    # Side-specific capsule groups
    for side in ["left", "right"]:
        capsule_term = f"{side} kidney capsule"
        capsule_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(capsule_term))
        capsule_exterior = fm.createFieldAnd(side_groups[side], is_exterior)
        capsule_group.getMeshGroup(mesh2d).addElementsConditional(capsule_exterior)
        side_kidney_groups[side][capsule_term] = capsule_group


def create_surface_and_edge_groups(fm, mesh2d, mesh1d, annotationGroups, region, side_groups, side_kidney_groups, is_exterior):
    """
    Create surface and edge annotation groups.
    """
    surface_types = ["anterior", "posterior", "lateral", "medial", "dorsal", "ventral", "juxtamedullary cortex"]
    surface_fields = {}

    # Create surface groups
    for surface_type in surface_types:
        if surface_type == "juxtamedullary cortex":
            cortex_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("cortex of kidney")).getGroup()
            medulla_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("renal medulla")).getGroup()
            surface_exterior = fm.createFieldAnd(medulla_group, cortex_group)
        else:
            base_group = getAnnotationGroupForTerm(annotationGroups, (surface_type, ""))
            surface_field = base_group.getGroup()
            surface_exterior = fm.createFieldAnd(surface_field, is_exterior)
            surface_fields[surface_type] = surface_field

        # General kidney surface group
        general_term = f"{surface_type} surface of kidney"
        general_surface_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(general_term))
        general_surface_group.getMeshGroup(mesh2d).addElementsConditional(surface_exterior)

        # Side-specific surface groups
        for side in ["left", "right"]:
            side_term = f"{surface_type} surface of {side} kidney"
            side_surface_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(side_term))
            side_surface_field = fm.createFieldAnd(surface_exterior, side_groups[side])
            side_surface_group.getMeshGroup(mesh2d).addElementsConditional(side_surface_field)
            side_kidney_groups[side][side_term] = side_surface_group

    # Create edge groups at dorsal-ventral intersection
    dorsal_ventral_border = fm.createFieldAnd(
        fm.createFieldAnd(surface_fields["dorsal"], surface_fields["ventral"]), is_exterior)

    edge_types = ["lateral", "medial"]

    for edge_type in edge_types:
        # General edge group
        edge_term = f"{edge_type} edge of kidney"
        edge_field = fm.createFieldAnd(surface_fields[edge_type], dorsal_ventral_border)
        general_edge_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(edge_term))
        general_edge_group.getMeshGroup(mesh1d).addElementsConditional(edge_field)

        # Side-specific edge groups
        for side in ["left", "right"]:
            side_edge_term = f"{edge_type} edge of {side} kidney"
            side_edge_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(side_edge_term))
            side_edge_field = fm.createFieldAnd(edge_field, side_groups[side])
            side_edge_group.getMeshGroup(mesh1d).addElementsConditional(side_edge_field)
            side_kidney_groups[side][side_edge_term] = side_edge_group


def remove_hidden_kidney_groups(show_kidneys, side_kidney_groups, annotationGroups):
    """
    Remove annotation groups for kidneys that should not be shown.
    """
    if not show_kidneys[0]:  # Left kidney
        for group in side_kidney_groups["left"].values():
            annotationGroups.remove(group)

    if not show_kidneys[1]:  # Right kidney
        for group in side_kidney_groups["right"].values():
            annotationGroups.remove(group)


def make_kidney_arcs(outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle, cortex_thickness,
                     arc_direction, number_of_elements):
    """
    Get kidney arcs on outer and sinus surface at a given arc angle up from +y toward +z
    These join in the form of a wave.
    :param outer_a: Kidney outer ellipsoid outer axis length in x direction.
    :param outer_b: Kidney outer ellipsoid outer axis length in y direction.
    :param outer_c: Kidney outer ellipsoid outer axis length in z direction.
    :param sinus_a: Kidney sinus ellipsoid outer axis length in x direction.
    :param sinus_b: Kidney sinus ellipsoid outer axis length in y direction.
    :param sinus_c: Kidney sinus ellipsoid outer axis length in z direction.
    :param cortex_thickness:
    :param hilum_angle: Angle from +x axis to hilum ellipse on sinus, PI/2 < hilum_angle < PI.
    :param arc_direction: Direction of arc plane, list [3] but expected in 2-3 plane.
    :param number_of_elements: Number of curve elements to sample around arc.
    :return: outer_x[], outer_d1[], sinus_x[], sinus_d1[]
    """
    assert (math.pi / 2.0) < hilum_angle < math.pi

    sinus_arc_angle = math.atan2(sinus_b * arc_direction[2], sinus_c * arc_direction[1])
    cos_sinus_arc_angle = math.cos(sinus_arc_angle)
    sin_sinus_arc_angle = math.sin(sinus_arc_angle)
    sinus_centre = [-sinus_a, 0.0, 0.0]
    sinus_axis1 = [sinus_a, 0.0, 0.0]
    sinus_axis2 = [0.0, sinus_b * cos_sinus_arc_angle, sinus_c * sin_sinus_arc_angle]
    sinus_x, sinus_d1 = sampleEllipsePoints(
        sinus_centre, sinus_axis1, sinus_axis2, 0.0, hilum_angle, number_of_elements)

    outer_arc_angle = math.atan2(outer_b * arc_direction[2], outer_c * arc_direction[1])
    cos_outer_arc_angle = math.cos(outer_arc_angle)
    sin_outer_arc_angle = math.sin(outer_arc_angle)
    outer_axis1 = [outer_a, 0.0, 0.0]
    outer_axis2 = [0.0, outer_b * cos_outer_arc_angle, outer_c * sin_outer_arc_angle]
    outer1_x, outer1_d1 = sampleEllipsePoints(
        [0.0, 0.0, 0.0], outer_axis1, outer_axis2, 0.0, math.pi / 2.0, number_of_elements)

    hilum_arc_angle = math.pi + math.acos(-sinus_x[-1][0] / outer_a)
    sin_hilum_arc_angle = math.sin(hilum_arc_angle)
    outer_return_b_mag = (
            (magnitude(outer_axis2) - magnitude([0.0, sinus_x[-1][1], sinus_x[-1][2]])) / (1.0 - sin_hilum_arc_angle))
    outer_return_axis2 = set_magnitude(outer_axis2, outer_return_b_mag)
    outer_return_centre = sub(outer_axis2, outer_return_axis2)

    outer2_x, outer2_d1 = sampleEllipsePoints(
        outer_return_centre, outer_axis1, outer_return_axis2, math.pi / 2.0, hilum_arc_angle, number_of_elements)

    outer_x, outer_d1 = sampleCubicHermiteCurves(
        outer1_x + outer2_x[1:], outer1_d1 + outer2_d1[1:], number_of_elements, arcLengthDerivatives=True)[0:2]

    return outer_x, outer_d1, sinus_x, sinus_d1


def build_kidney_capsule_octant(outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle,
                                element_counts, transition_count=1, nway_d_factor = 0.6):
    """

    :param outer_a:
    :param outer_b:
    :param outer_c:
    :param sinus_a:
    :param sinus_b:
    :param sinus_c:
    :param hilum_angle:
    :param element_counts:
    :param transition_count: Number of core transition elements
    :param nway_d_factor: Factor controlling derivative scale at n-way points.
    :return: HexTetrahedronMesh octant.
    """
    shell_count = 0
    around_count = 10
    around_node_count = around_count + 1
    over_count = 20
    over_node_count = over_count + 1
    delta_factor = 1.0E-5
    delta_factor2 = 2.0 * delta_factor


    outer_nx = []
    outer_nd1 = []
    outer_nd2 = []
    outer_nd12 = []
    sinus_nx = []
    sinus_nd1 = []
    sinus_nd2 = []
    sinus_nd12 = []

    arc_x, arc_d2 = sampleEllipsePoints([0.0, 0.0, 0.0], [0.0, outer_b, 0.0], [0.0, 0.0, outer_c], 0.0, math.pi / 2.0,
                                        around_count)

    for j in range(around_count + 1):
        outer_x, outer_d1, sinus_x, sinus_d1 = make_kidney_arcs(
            outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle, 0.0, arc_x[j], over_count)

        arc_dx = mult(arc_d2[j], delta_factor)
        outer_minus_x, outer_minus_d1, sinus_minus_x, sinus_minus_d1 = make_kidney_arcs(
            outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle, 0.0, sub(arc_x[j], arc_dx), over_count)
        outer_plus_x, outer_plus_d1, sinus_plus_x, sinus_plus_d1 = make_kidney_arcs(
            outer_a, outer_b, outer_c, sinus_a, sinus_b, sinus_c, hilum_angle, 0.0, add(arc_x[j], arc_dx), over_count)
        outer_d2 = [div(sub(outer_plus_x[i], outer_minus_x[i]), delta_factor2) for i in range(over_node_count)]
        outer_d12 = [div(sub(outer_plus_d1[i], outer_minus_d1[i]), delta_factor2) for i in range(over_node_count)]
        sinus_d2 = [div(sub(sinus_plus_x[i], sinus_minus_x[i]), delta_factor2) for i in range(over_node_count)]
        sinus_d12 = [div(sub(sinus_plus_d1[i], sinus_minus_d1[i]), delta_factor2) for i in range(over_node_count)]
        outer_nx += outer_x
        outer_nd1 += outer_d1
        outer_nd2 += outer_d2
        outer_nd12 += outer_d12
        sinus_nx += sinus_x
        sinus_nd1 += sinus_d1
        sinus_nd2 += sinus_d2
        sinus_nd12 += sinus_d12

    outer_surface = TrackSurface(over_count, around_count, outer_nx, outer_nd1, outer_nd2, outer_nd12)
    # node_identifier, element_identifier = outer_surface.generateMesh(region, node_identifier, element_identifier)

    sinus_surface = TrackSurface(over_count, around_count, sinus_nx, sinus_nd1, sinus_nd2, sinus_nd12)
    # node_identifier, element_identifier = sinus_surface.generateMesh(region, node_identifier, element_identifier)

    rim_count = shell_count + transition_count

    sample_curve_on_outer_surface = (
        lambda start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
               start_weight=None, end_weight=None, overweighting=1.0, end_transition=False:
        outer_surface.sampleCurve(
            start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
            start_weight, end_weight, overweighting, end_transition))
    move_x_to_outer_surface = lambda x: outer_surface.makeCoordinatesOnSurface(x)
    move_d_to_outer_surface = lambda x, d: outer_surface.makeCoordinatesAndDerivativeOnSurface(x, d)[1]

    sample_curve_on_sinus_surface = (
        lambda start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
               start_weight=None, end_weight=None, overweighting=1.0, end_transition=False:
        sinus_surface.sampleCurve(
            start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
            start_weight, end_weight, overweighting, end_transition))
    move_x_to_sinus_surface = lambda x: sinus_surface.makeCoordinatesOnSurface(x)
    move_d_to_sinus_surface = lambda x, d: sinus_surface.makeCoordinatesAndDerivativeOnSurface(x, d)[1]

    half_counts = [element_count // 2 for element_count in element_counts]
    diag_counts = [
        half_counts[0] + half_counts[1] - 2 * rim_count,
        half_counts[0] + half_counts[2] - 2 * rim_count,
        half_counts[1] + half_counts[2] - 2 * rim_count]
    octant = HexTetrahedronMesh(half_counts, diag_counts, nway_d_factor=nway_d_factor)

    box_counts = [half_counts[i] - rim_count for i in range(3)]

    abx, abd1, pe, pxi, psf = sampleCubicHermiteCurves(
        outer_nx[:over_node_count], outer_nd1[:over_node_count], diag_counts[0])
    abd2 = interpolateSampleCubicHermite(outer_nd2[:over_node_count], outer_nd12[:over_node_count], pe, pxi, psf)[0]

    acx, acd2, pe, pxi, psf = sampleCubicHermiteCurves(
        outer_nx[-over_node_count:], outer_nd1[-over_node_count:], diag_counts[1])
    acd1 = [[-d for d in d1] for d1 in interpolateSampleCubicHermite(
        outer_nd2[-over_node_count:], outer_nd12[-over_node_count:], pe, pxi, psf)[0]]

    bcx, bcd2, pe, pxi, psf = sampleCubicHermiteCurves(
        [outer_nx[over_count + i * over_node_count] for i in range(around_node_count)],
        [outer_nd2[over_count + i * over_node_count] for i in range(around_node_count)],
        diag_counts[2])
    bcd1 = interpolateSampleCubicHermite(
        [outer_nd1[over_count + i * over_node_count] for i in range(around_node_count)],
        [outer_nd12[over_count + i * over_node_count] for i in range(around_node_count)],
        pe, pxi, psf)[0]

    # set derivatives from other sampled sides of abc
    abd2[0] = acd2[0]
    abd2[-1] = bcd2[0]
    acd1[0] = abd1[0]
    acd1[-1] = [-d for d in bcd2[-1]]
    bcd1[0] = abd1[-1]
    bcd1[-1] = acd2[-1]

    triangle_abc = QuadTriangleMesh(
        box_counts[0], box_counts[1], box_counts[2],
        sample_curve_on_outer_surface, move_x_to_outer_surface, move_d_to_outer_surface,
        nway_d_factor=nway_d_factor)
    triangle_abc.set_edge_parameters12(abx, abd1, abd2)
    triangle_abc.set_edge_parameters13(acx, acd1, acd2)
    triangle_abc.set_edge_parameters23(bcx, bcd1, bcd2)
    triangle_abc.build()
    # GRC review:
    normal_size = 0.1
    triangle_abc.assign_d3(lambda tx, td1, td2: set_magnitude(cross(td1, td2), normal_size))
    octant.set_triangle_abc(triangle_abc)

    # extract actual derivatives calculated on edges of triangle abc
    abx, abd1, abd2, abd3 = triangle_abc.get_edge_parameters12()
    acx, acd1, acd2, acd3 = triangle_abc.get_edge_parameters13()
    bcx, bcd1, bcd2, bcd3 = triangle_abc.get_edge_parameters23()

    # get parameters from a, b, c back to origin
    size_ao = outer_a / half_counts[0]
    aox = [[size_ao * i, 0.0, 0.0] for i in range(half_counts[0], -1, -1)]
    aod2 = [[-size_ao, 0.0, 0.0]] * (half_counts[0] + 1)
    aod1 = [[0.0, normal_size, 0.0]] * (half_counts[0] + 1)
    aod3 = [[0.0, 0.0, normal_size]] * (half_counts[0] + 1)

    obx, obd2, pe, pxi, psf = sampleCubicHermiteCurves(
        sinus_nx[:over_node_count], sinus_nd1[:over_node_count], half_counts[1])
    obd3 = [[0.0, 0.0, normal_size]] * (half_counts[1] + 1)
    obd1 = [set_magnitude([-d2[1], d2[0], 0.0], normal_size) for d2 in obd2]
    box = list(reversed(obx))
    bod1 = list(reversed(obd1))
    bod2 = [[-d for d in d2] for d2 in reversed(obd2)]
    bod3 = list(reversed(obd3))
    del obx, obd1, obd2, obd3

    ocx, ocd2, pe, pxi, psf = sampleCubicHermiteCurves(
        sinus_nx[-over_node_count:], sinus_nd1[-over_node_count:], half_counts[2])
    ocd3 = [[0.0, -normal_size, 0.0]] * (half_counts[2] + 1)
    ocd1 = [set_magnitude([-d2[2], 0.0, d2[0]], normal_size) for d2 in ocd2]
    cox = list(reversed(ocx))
    cod1 = list(reversed(ocd1))
    cod2 = [[-d for d in d2] for d2 in reversed(ocd2)]
    cod3 = list(reversed(ocd3))
    del ocx, ocd1, ocd2, ocd3

    # fix known derivatives including zero derivative at pole
    aod1[0] = abd1[0]
    aod1[-1] = [-d for d in bod2[-1]]
    aod3[0] = abd2[0]
    aod3[-1] = bod3[-1]
    bod1[0] = abd1[-1]
    bod1[-1] = aod2[-1]
    bod3[0] = bcd2[0]
    bod3[-1] = [-d for d in cod2[-1]]
    cod1[0] = acd2[-1]
    cod1[-1] = aod2[-1]
    cod3[0] = bcd2[-1]
    cod3[-1] = bod2[-1]

    # make inner surface triangle abo
    triangle_abo = QuadTriangleMesh(
        box_counts[0], box_counts[1], rim_count, sampleHermiteCurve, nway_d_factor=nway_d_factor)
    mabd3 = [[-d for d in d3] for d3 in abd3]
    triangle_abo.set_edge_parameters12(abx, abd1, mabd3, abd2)
    triangle_abo.set_edge_parameters13(aox, aod1, aod2, aod3)
    triangle_abo.set_edge_parameters23(box, bod1, bod2, bod3)
    triangle_abo.build()
    # GRC review:
    triangle_abo.assign_d3(lambda tx, td1, td2: set_magnitude(cross(td1, td2), 0.1))
    octant.set_triangle_abo(triangle_abo)
    bod1 = triangle_abo.get_edge_parameters23()[1]

    # make inner surface triangle aco
    triangle_aco = QuadTriangleMesh(
        box_counts[0], box_counts[2], rim_count, sampleHermiteCurve, nway_d_factor=nway_d_factor)
    acmd3 = [[-d for d in d3] for d3 in acd3]
    acmd1 = [[-d for d in d1] for d1 in acd1]
    aomd1 = [[-d for d in d1] for d1 in aod1]
    triangle_aco.set_edge_parameters12(acx, acd2, acmd3, acmd1)
    triangle_aco.set_edge_parameters13(aox, aod3, aod2, aomd1)
    triangle_aco.set_edge_parameters23(cox, cod1, cod2, cod3)
    triangle_aco.build()
    # GRC review:
    triangle_aco.assign_d3(lambda tx, td1, td2: set_magnitude(cross(td1, td2), normal_size))
    octant.set_triangle_aco(triangle_aco)
    cod1 = triangle_aco.get_edge_parameters23()[1]

    # make inner surface bco
    triangle_bco = QuadTriangleMesh(
        box_counts[1], box_counts[2], rim_count,
        sample_curve_on_sinus_surface, move_x_to_sinus_surface, move_d_to_sinus_surface,
        nway_d_factor=nway_d_factor)
    # GRC needed?
    _, bcd1, _, bcd3 = triangle_abc.get_edge_parameters23()
    # substitute end derivatives for which magnitude is known:
    # bcd3[0] = bod3[0]
    # bcd3[-1] = cod2[-1]
    bcmd1 = [[-d for d in d1] for d1 in bcd1]
    bcmd3 = [[-d for d in d3] for d3 in bcd3]
    triangle_bco.set_edge_parameters12(bcx, bcd2, bcmd3, bcmd1)
    bomd1 = [[-d for d in d1] for d1 in bod1]
    # bomd3 = [[-d for d in d3] for d3 in bod3]
    # triangle_bco.set_edge_parameters13(box, bod2, bomd1, bomd3)
    triangle_bco.set_edge_parameters13(box, bod3, bod2, bomd1)
    comd1 = [[-d for d in d1] for d1 in cod1]
    triangle_bco.set_edge_parameters23(cox, cod3, cod2, comd1)
    triangle_bco.build()
    # GRC review:
    triangle_bco.assign_d3(lambda tx, td1, td2: set_magnitude(cross(td1, td2), normal_size))
    octant.set_triangle_bco(triangle_bco)

    # _, aco_cod1, _, aco_cod3 = triangle_aco.get_edge_parameters23()

    # # fix known derivatives
    # aod2[0] = abd1[0]
    # aod3[0] = [-d for d in abd2[0]]
    #
    # aox, aod2, aod1 = sampleHermiteCurve(
    #     ext_axis1, ext_axis_md1, abd1[0], ext_origin, ext_axis_md1, axis_d2, elements_count=half_counts[0])
    # box, bod2, bod3 = sampleHermiteCurve(
    #     axis2, axis_md2, bcd2[0], origin, axis_md2, axis_d3, elements_count=half_counts[1])
    # bod1 = [abd1[-1]] + [axis_md1] * (len(box) - 1)
    # cox, cod2, cod1 = sampleHermiteCurve(
    #     ext_axis3, axis_md3, [-d for d in acd1[-1]], ext_origin, axis_md3, bod2[-1],
    #     elements_count=half_counts[2])

    #
    # aox, aod1, aod2, aod3 = crx[l], crd1[l], crd2[l], crd3[l]
    # box, bod1, bod2, bod3 = crx[lp], crd1[lp], crd2[lp], crd3[lp]
    # # get last tube row core centre coordinates
    # n1 = self._elementsCountCoreBoxMajor // 2
    # n2 = dome_ix
    # n3 = self._elementsCountCoreBoxMinor // 2
    # ccx = self._boxCoordinates[0][n2][n1][n3]
    # ccd1 = self._boxCoordinates[1][n2][n1][n3]
    # ccd2 = self._boxCoordinates[2][n2][n1][n3]
    # ccd3 = self._boxCoordinates[3][n2][n1][n3]
    # # get dome pole coordinates. Note d1 is zero at the pole so use d2 in separate quadrants
    # core_p = 1 if self._shell_count else 0
    # dpx = plx[core_p][0][-1]
    # dpd1 = pld2[core_p][0][-1]
    # dpd2 = self._pathParametersList[core_p][1][dome_ix]  # direction only
    # dpd3 = pld2[core_p][1][-1]
    # if dome0:
    #     ccd2 = [-d for d in ccd2]
    #     dpd2 = [-d for d in dpd2]
    # else:
    #     ccd3 = [-d for d in ccd3]
    # # sample inner line from dome pole back to last tube line core centre
    # bx = ccx
    # bd = mult(ccd2, -half_counts[2])
    # ax = dpx
    # ad = computeCubicHermiteStartDerivative(ax, [-d for d in dpd2], bx, bd)
    # cox = []
    # cod1 = []
    # cod2 = []
    # cod3 = []
    # for n in range(half_counts[2] + 1):
    #     xi = n / half_counts[2]
    #     # future: use cubic expression for twist in d1, d3
    #     cox.append(interpolateCubicHermite(ax, ad, bx, bd, xi))
    #     cod1.append(linearlyInterpolateVectors(dpd1, ccd1, xi))
    #     cod2.append(div(interpolateCubicHermiteDerivative(ax, ad, bx, bd, xi), half_counts[2]))
    #     cod3.append(linearlyInterpolateVectors(dpd3, ccd3, xi))
    #
    # # make inner surface triangle 1-2-origin
    # triangle_abo = QuadTriangleMesh(
    #     box_counts[0], box_counts[1], rim_count, sampleHermiteCurve, nway_d_factor=self._nway_d_factor)
    # abd3 = [[-d for d in d3] for d3 in abd3]
    # triangle_abo.set_edge_parameters12(abx, abd1, abd3, abd2)
    # aomd1 = [[-d for d in d1] for d1 in aod1]
    # triangle_abo.set_edge_parameters13(aox, aomd1, aod3, aod2)
    # bomd1 = [[-d for d in d1] for d1 in bod1]
    # triangle_abo.set_edge_parameters23(box, bomd1, bod3, bod2)
    # triangle_abo.build()
    # triangle_abo.assign_d3(lambda tx, td1, td2: normalize(cross(td1, td2)))
    # core_octant.set_triangle_abo(triangle_abo)
    # # extract actual derivatives calculated on edges of triangle abo
    # _, _, aod3, _ = triangle_abo.get_edge_parameters13()
    # _, bomd1, bod3, _ = triangle_abo.get_edge_parameters23()
    # bod1 = [[-d for d in d1] for d1 in bomd1]
    #
    # # make inner surface triangle 1-3-origin
    # triangle_aco = QuadTriangleMesh(
    #     box_counts[0], box_counts[2], rim_count, sampleHermiteCurve, nway_d_factor=self._nway_d_factor)
    # acd3 = [[-d for d in d3] for d3 in triangle_abc.get_edge_parameters13()[3]]
    # acmd1 = [[-d for d in d1] for d1 in acd1]
    # triangle_aco.set_edge_parameters12(acx, acd2, acd3, acmd1)
    # triangle_aco.set_edge_parameters13(aox, aod2, aod3, aod1)
    # comd1 = [[-d for d in d1] for d1 in cod1]
    # comd3 = [[-d for d in d3] for d3 in cod3]
    # if l == 0:
    #     use_cod1 = cod1
    #     use_cod3 = cod3
    # elif l == 1:
    #     use_cod1 = cod3
    #     use_cod3 = comd1
    # elif l == 2:
    #     use_cod1 = comd1
    #     use_cod3 = comd3
    # else:  # if l == 3:
    #     use_cod1 = comd3
    #     use_cod3 = cod1
    # triangle_aco.set_edge_parameters23(cox, use_cod1, cod2, use_cod3)
    #
    # triangle_aco.build()
    # triangle_aco.assign_d3(lambda tx, td1, td2: normalize(cross(td1, td2)))
    # core_octant.set_triangle_aco(triangle_aco)
    # _, aco_cod1, _, aco_cod3 = triangle_aco.get_edge_parameters23()
    #
    # # make inner surface 2-3-origin
    # triangle_bco = QuadTriangleMesh(
    #     box_counts[1], box_counts[2], rim_count, sampleHermiteCurve, nway_d_factor=self._nway_d_factor)
    # _, bcd1, _, bcd3 = triangle_abc.get_edge_parameters23()
    # # substitute end derivatives for which magnitude is known:
    # # bcd3[0] = bod3[0]
    # # bcd3[-1] = cod2[-1]
    # bcmd1 = [[-d for d in d1] for d1 in bcd1]
    # bcmd3 = [[-d for d in d3] for d3 in bcd3]
    # triangle_bco.set_edge_parameters12(bcx, bcd2, bcmd3, bcmd1)
    # triangle_bco.set_edge_parameters13(box, bod2, bod3, bod1)
    # triangle_bco.set_edge_parameters23(cox, aco_cod3, cod2, [[-d for d in d3] for d3 in aco_cod1])
    # triangle_bco.build()
    # triangle_bco.assign_d3(lambda tx, td1, td2: normalize(cross(td1, td2)))
    # core_octant.set_triangle_bco(triangle_bco)
    #

    octant.build_interior()
    return octant
