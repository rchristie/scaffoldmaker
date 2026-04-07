"""
Utility class for defining network meshes from 1-D connectivity and lateral axes, with continuity control.
"""
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.field import Field
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup


class MeshGenerateData:
    """
    Data for passing to generate mesh functions.
    Maintains Zinc region, field, node and element numbering, node layout map, and output annotation groups.
    Derive from this class to pass additional data.
    """

    def __init__(self, region, meshDimension, coordinateFieldName="coordinates",
                 startNodeIdentifier=1, startElementIdentifier=1):
        self._region = region
        self._fieldmodule = region.getFieldmodule()
        self._fieldcache = self._fieldmodule.createFieldcache()
        self._meshDimension = meshDimension
        self._mesh = self._fieldmodule.findMeshByDimension(meshDimension)
        self._nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._coordinates = find_or_create_field_coordinates(self._fieldmodule, coordinateFieldName)
        self._nodeIdentifier = startNodeIdentifier
        self._elementIdentifier = startElementIdentifier
        self._annotationGroups = []  # list of AnnotationGroup to return for mesh's scaffold
        self._annotationGroupMap = {}  # map from annotation term (name, ontId) to AnnotationGroup in output region

    def getCoordinates(self):
        """
        :return: Zinc Finite Element coordinate field being defined.
        """
        return self._coordinates

    def getFieldcache(self):
        """
        :return: Zinc Fieldcache for assigning field parameters with.
        """
        return self._fieldcache

    def getFieldmodule(self):
        """
        :return: Zinc Fieldmodule being generated in.
        """
        return self._fieldmodule

    def getMesh(self):
        """
        :return: Zinc Mesh for elements being built.
        """
        return self._mesh

    def getMeshDimension(self):
        """
        :return: Dimension of elements being built.
        """
        return self._meshDimension

    def getNodes(self):
        """
        :return: Zinc Nodeset for nodes being built.
        """
        return self._nodes

    def getRegion(self):
        return self._region

    def getNodeElementIdentifiers(self):
        """
        Get next node and element identifiers without incrementing, to call at end of generation.
        :return: Next node identifier, next element identifier.
        """
        return self._nodeIdentifier, self._elementIdentifier

    def setNodeElementIdentifiers(self, nodeIdentifier, elementIdentifier):
        """
        Set next node and element identifiers after generating objects with external code.
        """
        self._nodeIdentifier = nodeIdentifier
        self._elementIdentifier = elementIdentifier

    def nextNodeIdentifier(self):
        nodeIdentifier = self._nodeIdentifier
        self._nodeIdentifier += 1
        return nodeIdentifier

    def nextElementIdentifier(self):
        elementIdentifier = self._elementIdentifier
        self._elementIdentifier += 1
        return elementIdentifier

    def getRegion(self):
        return self._region

    def getAnnotationGroups(self):
        return self._annotationGroups

    def getOrCreateAnnotationGroup(self, annotationTerm):
        annotationGroup = self._annotationGroupMap.get(annotationTerm)
        if not annotationGroup:
            annotationGroup = AnnotationGroup(self._region, annotationTerm)
            self._annotationGroups.append(annotationGroup)
            self._annotationGroupMap[annotationTerm] = annotationGroup
        return annotationGroup

    def _getAnnotationMeshGroup(self, annotationTerm):
        """
        Get mesh group to add elements to for term.
        :param annotationTerm: Annotation term (name, ontId).
        :return: Zinc MeshGroup.
        """
        annotationGroup = self.getOrCreateAnnotationGroup(annotationTerm)
        return annotationGroup.getMeshGroup(self._mesh)

    def getAnnotationMeshGroups(self, annotationTerms):
        """
        Get mesh groups for all annotation terms to add segment elements to, creating as needed.
        :param annotationTerms: List of annotation terms (name, ontId).
        :return: List of Zinc MeshGroup.
        """
        return [self._getAnnotationMeshGroup(annotationTerm) for annotationTerm in annotationTerms]
