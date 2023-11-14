"""
Constructs a 1-D network layout mesh with specifiable structure.
"""
from cmlibs.maths.vectorops import cross, magnitude, mult, normalize, sub
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.scene import scene_get_selection_group
from cmlibs.zinc.field import Field, FieldGroup
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.interpolation import smoothCubicHermiteCrossDerivativesLine
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.zinc_utils import clearRegion, get_nodeset_field_parameters, \
    get_nodeset_path_field_parameters, make_nodeset_derivatives_orthogonal, scale_nodeset_derivatives, \
    set_nodeset_field_parameters, setPathParameters
from enum import Enum
import math


class MeshType_1d_network_layout1(Scaffold_base):
    """
    Defines branching network layout with side dimensions.
    """

    parameterSetStructureStrings = {
        "Default": "1-2",
        "Bifurcation": "1-2,2-3,2.2-4",
        "Converging bifurcation": "1-3.1,2-3.2,3.3-4",
        "Sphere cube": "1.1-2.1,1.2-3.1,1.3-4.1,2.2-5.2,2.3-6.1,3.2-6.2,3.3-7.1,4.2-7.2,4.3-5.1,5.3-8.1,6.3-8.2,7.3-8.3"
    }

    @classmethod
    def getName(cls):
        return "1D Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default", "Bifurcation", "Converging bifurcation", "Sphere cube"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = parameterSetName
        options["Structure"] = cls.parameterSetStructureStrings[parameterSetName]
        options["Define inner coordinates"] = False  # can be overridden by parent scaffold
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            #  "Base parameter set"  # Hidden.
            #  "Structure"  # Hidden so must edit via interactive function.
            #  "Define inner coordinates"  # Hidden as enabled by parent scaffold.
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, NetworkMesh
        """
        parameterSetName = options['Base parameter set']
        structure = options["Structure"]
        defineInnerCoordinates = options["Define inner coordinates"]
        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule).castFiniteElement()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        fieldcache = fieldmodule.createFieldcache()
        if "Sphere cube" in parameterSetName:
            # edit node parameters
            sphereRadius = 0.5
            tubeRadius = 0.1
            edgeAngle = 2.0 * math.asin(math.sqrt(1.0 / 3.0))
            # get x and d3
            cx = []
            cd3 = []
            for i in range(4):
                angleUp = [0.0, edgeAngle, math.pi - edgeAngle, math.pi][i]
                cosAngleUp = math.cos(angleUp)
                sinAngleUp = math.sin(angleUp)
                z = -sphereRadius * cosAngleUp
                zRadius = sphereRadius * sinAngleUp
                jLimit = 1 if i in [0, 3] else 3
                for j in range(jLimit):
                    angleAround = math.radians(120.0 * ((j - 0.5) if (i == 2) else j))
                    cosAngleAround = math.cos(angleAround)
                    sinAngleAround = math.sin(angleAround)
                    px = [zRadius * cosAngleAround, zRadius * sinAngleAround, z]
                    cx.append(px)
                    cd3.append(mult(normalize(px), tubeRadius))
            # get d1, d2, d13
            cd1 = []
            cd2 = []
            cd13 = []
            for n in range(8):
                cd1.append([])
                cd2.append([])
                cd13.append([])
            edgeArcLength = sphereRadius * edgeAngle
            for networkSegment in networkMesh.getNetworkSegments():
                networkNodes = networkSegment.getNetworkNodes()
                nodeIndexes = [networkNode.getNodeIdentifier() - 1 for networkNode in networkNodes]
                delta = sub(cx[nodeIndexes[1]], cx[nodeIndexes[0]])
                for ln in range(2):
                    d3 = cd3[nodeIndexes[ln]]
                    d2 = mult(normalize(cross(d3, delta)), tubeRadius)
                    d1 = mult(normalize(cross(d2, d3)), edgeArcLength)
                    cd1[nodeIndexes[ln]].append(d1)
                    cd2[nodeIndexes[ln]].append(d2)
                    cd13[nodeIndexes[ln]].append(mult(d1, tubeRadius))
            # fix the one node out of order:
            for d in [cd1[4], cd2[4]]:
                d[0:2] = [d[1], d[0]]
            nodeId = 1
            for n in range(8):
                node = nodes.findNodeByIdentifier(nodeId)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
                for v in range(3):
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, v + 1, cd1[n][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, v + 1, cd2[n][v])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, v + 1, cd3[n])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, v + 1, cd13[n][v])
                nodeId += 1

            if defineInnerCoordinates:
                # copy coordinates to inner coordinates via in-memory model file
                coordinates.setName("inner coordinates")  # temporarily rename
                sir = region.createStreaminformationRegion()
                srm = sir.createStreamresourceMemory()
                region.write(sir)
                result, buffer = srm.getBuffer()
                coordinates.setName("coordinates")  # restore name before reading inner coordinates back in
                sir = region.createStreaminformationRegion()
                srm = sir.createStreamresourceMemoryBuffer(buffer)
                region.read(sir)
                functionOptions = {
                    "Mode": {"Wall proportion": True, "Wall thickness": False},
                    "Value": 0.5}
                cls.assignInnerCoordinates(region, options, networkMesh, functionOptions, editGroupName=None)

        return [], networkMesh

    @classmethod
    def editStructure(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Edit structure safely, to prevent accidental changes.
        Copies functionOptions["Structure"] to options["Structure"] and regenerates with
        default geometric coordinates.
        :param region: Region containing model to clear and re-generate.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Contents replaced.
        :param functionOptions: functionOptions["Structure"] contains new structure string.
        :param editGroupName: Name of Zinc group to put edited nodes in. Cleared.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        with ChangeManager(fieldmodule):
            clearRegion(region)
            structure = options["Structure"] = functionOptions["Structure"]
            networkMesh.build(structure)
            networkMesh.create1DLayoutMesh(region)
        return True, False  # settings changed, nodes not changed (since reset to original coordinates)

    class InnerCoordinatesMode(Enum):
        """
        Controls how inner
        """
        PROPORTIONAL_WALL_THICKNESS = 1  # same
        ABSOLUTE_WALL_THICKNESS = 2

    @classmethod
    def assignInnerCoordinates(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Assign inner coordinates from outer coordinates with wall thickness either a proportion of the
        outer side derivatives, or with an absolute wall thickness.
        If elements are selected, applies
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Unused.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        innerCoordinates = fieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if not innerCoordinates.isValid():
            print("Assign inner coordinates:  No inner coordinates defined")
            return None, None
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        assignMode = None
        if functionOptions["Mode"]["Wall proportion"]:
            assignMode = cls.InnerCoordinatesMode.PROPORTIONAL_WALL_THICKNESS
        elif functionOptions["Mode"]["Wall thickness"]:
            assignMode = cls.InnerCoordinatesMode.ABSOLUTE_WALL_THICKNESS
        else:
            assert False, "assignInnerCoordinates. Invalid mode"
        wallThicknessValue = functionOptions["Value"]
        selectionGroup = scene_get_selection_group(region.getScene(), inherit_root_region=region.getRoot())
        selectionMeshGroup = None
        mesh1d = fieldmodule.findMeshByDimension(1)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        if selectionGroup:
            selectionMeshGroup = selectionGroup.getMeshGroup(mesh1d)
            if not selectionMeshGroup.isValid():
                print("Assign inner coordinates:  Selection contains no elements. Clear it to assign globally.")
                return None, None
        valueLabels = [
            Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        with ChangeManager(fieldmodule):
            # get all node parameters (from selection if any)
            useNodeset = nodeset
            if selectionGroup:
                tmpGroup = fieldmodule.createFieldGroup()
                tmpGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
                tmpMeshGroup = tmpGroup.createMeshGroup(mesh1d)
                tmpMeshGroup.addElementsConditional(selectionGroup)
                useNodeset = tmpGroup.getNodesetGroup(nodes)
                print("selection nodeset size", useNodeset.getSize())
                del tmpMeshGroup
                del tmpGroup
            _, nodeParameters = get_nodeset_field_parameters(useNodeset, coordinates, valueLabels)

            nodeIdentifierIndexes = {}
            modifyVersions = []
            for n in range(len(nodeParameters)):
                nodeIdentifierIndexes[nodeParameters[n][0]] = n
                versionsCount = len(nodeParameters[n][1][1])
                modifyVersions.append([False] * versionsCount)

            networkSegments = networkMesh.getNetworkSegments()
            for networkSegment in networkSegments:
                nodeIdentifiers = networkSegment.getNodeIdentifiers()
                nodeVersions = networkSegment.getNodeVersions()
                for n in range(len(nodeIdentifiers)):
                    nodeIndex = nodeIdentifierIndexes.get(nodeIdentifiers[n])
                    print("Node identifier", nodeIdentifiers[n], "index", nodeIndex, "version", nodeVersions[n])
                    if nodeIndex is not None:
                        modifyVersions[nodeIndex][nodeVersions[n] - 1] = True

            proportion = 1.0 - wallThicknessValue
            for n in range(len(nodeParameters)):
                modifyVersion = modifyVersions[n]
                versionsCount = len(modifyVersion)
                nNodeParameters = nodeParameters[n][1]
                for v in range(versionsCount):
                    if modifyVersion[v]:
                        if assignMode == cls.InnerCoordinatesMode.PROPORTIONAL_WALL_THICKNESS:
                            for d in range(2, 6):
                                nNodeParameters[d][v] = mult(nNodeParameters[d][v], proportion)
                        elif assignMode == cls.InnerCoordinatesMode.ABSOLUTE_WALL_THICKNESS:
                            for dd in range(2):
                                mag = magnitude(nNodeParameters[2 + 2 * dd][v])
                                if abs(mag) > 0.0:
                                    proportion = (mag - wallThicknessValue) / mag
                                    for d in [2 + 2 * dd, 3 + 2 * dd]:
                                        nNodeParameters[d][v] = mult(nNodeParameters[d][v], proportion)

            set_nodeset_field_parameters(useNodeset, innerCoordinates, valueLabels, nodeParameters, editGroupName)
            del useNodeset

        return False, True  # settings not changed, nodes changed

    @classmethod
    def makeSideDerivativesNormal(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Make side directions normal to d1 and each other. Works for all versions.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Unused.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        makeD2Normal = functionOptions['Make D2 normal']
        makeD3Normal = functionOptions['Make D3 normal']
        if not (makeD2Normal or makeD3Normal):
            return False, False
        make_nodeset_derivatives_orthogonal(nodeset, coordinates, makeD2Normal, makeD3Normal, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def smoothSideCrossDerivatives(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Smooth side cross derivatives giving rate of change of side directions d2, d3 w.r.t. d1.
        Note: only works for a single path with version 1.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from.
        Used to determine connected paths for smoothing.
        :param functionOptions: Which side derivatives to smooth.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        smoothD12 = functionOptions["Smooth D12"]
        smoothD13 = functionOptions["Smooth D13"]
        if not (smoothD12 or smoothD13):
            return False, False
        valueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1]
        if smoothD12:
            valueLabels += [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2]
        if smoothD13:
            valueLabels += [Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]
        fieldmodule = region.getFieldmodule()
        parameters = get_nodeset_path_field_parameters(
            fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES),
            fieldmodule.findFieldByName('coordinates'),
            valueLabels)
        x = parameters[0]
        d1 = parameters[1]
        modifyParameters = []
        modifyValueLabels = []
        if smoothD12:
            d12 = smoothCubicHermiteCrossDerivativesLine(x, d1, parameters[2], parameters[3])
            modifyParameters.append(d12)
            modifyValueLabels.append(Node.VALUE_LABEL_D2_DS1DS2)
        if smoothD13:
            d13 = smoothCubicHermiteCrossDerivativesLine(x, d1, parameters[-2], parameters[-1])
            modifyParameters.append(d13)
            modifyValueLabels.append(Node.VALUE_LABEL_D2_DS1DS3)
        setPathParameters(region, modifyValueLabels, modifyParameters, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def scaleSideDerivatives(cls, region, options, networkMesh, functionOptions, editGroupName):
        """
        Make side directions normal to d1 and each other. Works for all versions.
        :param region: Region containing model to change parameters of.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param networkMesh: The NetworkMesh construction object model was created from. Unused.
        :param functionOptions: Which side directions to make normal.
        :param editGroupName: Name of Zinc group to put edited nodes in.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
        nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        d2Scale = functionOptions['D2 scale']
        d3Scale = functionOptions['D3 scale']
        if (d2Scale == 0.0) or (d3Scale == 0.0):
            return False, False
        valueLabels = []
        scalings = []
        if d2Scale != 1.0:
            valueLabels += [Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2]
            scalings += [d2Scale, d2Scale]
        if d3Scale != 1.0:
            valueLabels += [Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]
            scalings += [d3Scale, d3Scale]
        if not scalings:
            return False, False
        scale_nodeset_derivatives(nodeset, coordinates, valueLabels, scalings, editGroupName)
        return False, True  # settings not changed, nodes changed

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Supply client with functions for smoothing path parameters.
        """
        return Scaffold_base.getInteractiveFunctions() + [
            ("Edit structure...",
                {"Structure": None},  # None = take value from options
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.editStructure(region, options, networkMesh, functionOptions, editGroupName)),
            ("Assign inner coordinates...",
                {"Mode": {"Wall proportion": True, "Wall thickness": False},
                 "Value": 0.5},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.assignInnerCoordinates(region, options, networkMesh, functionOptions, editGroupName)),
            ("Make side derivatives normal...",
                {"Make D2 normal": True,
                 "Make D3 normal": True},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.makeSideDerivativesNormal(region, options, networkMesh, functionOptions, editGroupName)),
            ("Smooth side cross derivatives...",
                {"Smooth D12": True,
                 "Smooth D13": True},
                lambda region, options, networkMesh, functionOptions, editGroupName:
                    cls.smoothSideCrossDerivatives(region, options, networkMesh, functionOptions, editGroupName)),
            ("Scale side derivatives...",
             {"D2 scale": 1.0,
              "D3 scale": 1.0},
             lambda region, options, networkMesh, functionOptions, editGroupName:
                cls.scaleSideDerivatives(region, options, networkMesh, functionOptions, editGroupName)),
        ]
