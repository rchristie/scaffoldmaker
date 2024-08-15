"""
Constructs a 1-D linear graph.
"""
from cmlibs.maths.vectorops import add, mult
from cmlibs.utils.zinc.field import find_or_create_field_coordinates, find_or_create_field_finite_element
from cmlibs.utils.zinc.finiteelement import get_maximum_element_identifier, get_maximum_node_identifier
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_1d_linear_graph1(Scaffold_base):
    """
    Defines 1-D linear graph with radius, suitable for segmentation mock-up.
    """

    @classmethod
    def getName(cls):
        return "1D Linear Graph 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Mesh buffer"] = ""
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            # "Mesh"  created interactively
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
        # nothing to do: entirely created interactively
        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        meshBuffer = options["Mesh buffer"]
        if meshBuffer:
            sir = region.createStreaminformationRegion()
            srm = sir.createStreamresourceMemoryBuffer(meshBuffer)
            region.read(sir)
        return [], None

    @classmethod
    def addSegment(cls, region, options, constructionObject, functionOptions, editGroupName):
        """
        Add a sequence of nodes and elements from start to end coordinates with radius.
        If start node or end node supplied, coordinates are taken from it.
        :param region: Region containing model to clear and re-generate.
        :param options: The scaffold settings used to create the original model, pre-edits.
        :param constructionObject: The construction object model was created from. Always None.
        :param functionOptions: functionOptions["Structure"] contains new structure string.
        :param editGroupName: Name of Zinc group to put edited nodes in. Cleared.
        :return: boolean indicating if settings changed, boolean indicating if node parameters changed.
        """
        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES);
        mesh1d = fieldmodule.findMeshByDimension(1)
        nodesCount = functionOptions["Number of nodes"]
        with ChangeManager(fieldmodule):
            coordinates = find_or_create_field_coordinates(fieldmodule)
            radius = find_or_create_field_finite_element(fieldmodule, "radius", 1)
            fieldcache = fieldmodule.createFieldcache()
            startNode = nodes.findNodeByIdentifier(functionOptions["Start node"])
            if startNode.isValid():
                fieldcache.setNode(startNode)
                result, startCoordinates = coordinates.evaluateReal(fieldcache, 3)
            else:
                startNode = None
                startCoordinates = functionOptions["Start x, y, z"]
                if len(startCoordinates) > 3:
                    startCoordinates = startCoordinates[:3]
                else:
                    while len(startCoordinates) < 3:
                        startCoordinates.append(0.0)
            endNode = nodes.findNodeByIdentifier(functionOptions["End node"])
            if endNode.isValid():
                fieldcache.setNode(endNode)
                result, endCoordinates = coordinates.evaluateReal(fieldcache, 3)
            else:
                endNode = None
                endCoordinates = functionOptions["End x, y, z"]
                if len(endCoordinates) > 3:
                    endCoordinates = startCoordinates[:3]
                else:
                    while len(endCoordinates) < 3:
                        endCoordinates.append(0.0)
            radiusValue = functionOptions["Radius"]
            nodeIdentifier = max(1, get_maximum_node_identifier(nodes) + 1)
            elementIdentifier = max(1, get_maximum_element_identifier(mesh1d) + 1)
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            if radius:
                 nodetemplate.defineField(radius)
            elementtemplate = mesh1d.createElementtemplate()
            elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
            linearbasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
            eft = mesh1d.createElementfieldtemplate(linearbasis)
            elementtemplate.defineField(coordinates, -1, eft)
            if radius:
                 elementtemplate.defineField(radius, -1, eft)
            lastNode = None

            xi_scale = 1.0 / max(1, nodesCount - 1)
            for n in range(nodesCount):
                xi = n * xi_scale
                if startNode and (n == 0):
                    node = startNode
                elif endNode and (n == (nodesCount - 1)):
                    node = endNode
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    x = add(mult(startCoordinates, 1.0 - xi), mult(endCoordinates, xi))
                    fieldcache.setNode(node)
                    coordinates.assignReal(fieldcache, x)
                    radius.assignReal(fieldcache, radiusValue)
                    nodeIdentifier += 1
                if lastNode:
                    element = mesh1d.createElement(elementIdentifier, elementtemplate)
                    element.setNode(eft, 1, lastNode)
                    element.setNode(eft, 2, node)
                    elementIdentifier += 1
                lastNode = node
        editGroup = fieldmodule.findFieldByName(editGroupName).castGroup()
        if editGroup.isValid():
            editGroup.clear()
            editGroup.setManaged(False)
        del editGroup
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        region.write(sir)
        result, meshBuffer = srm.getBuffer()
        options["Mesh buffer"] = meshBuffer
        return True, True  # settings changed, nodes changed

    @classmethod
    def getInteractiveFunctions(cls):
        """
        Supply client with functions for smoothing path parameters.
        """
        return Scaffold_base.getInteractiveFunctions() + [
            ("Add segment...", {
                "Number of nodes": 2,
                "Start node": 0,
                "Start x, y, z": [0.0, 0.0, 0.0],
                "End node": 0,
                "End x, y, z": [1.0, 0.0, 0.0],
                "Radius": 0.1},
                lambda region, options, constructionObject, functionOptions, editGroupName:
                    cls.addSegment(region, options, constructionObject, functionOptions, editGroupName))
        ]
