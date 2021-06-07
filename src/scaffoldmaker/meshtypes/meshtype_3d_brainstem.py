"""
Brainstem mesh using a tapered cylinder
"""

from __future__ import division
import copy
from opencmiss.zinc.element import Element
from opencmiss.zinc.node import Node
from opencmiss.zinc.field import Field, FieldFindMeshLocation
from opencmiss.utils.zinc.field import Field, findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldNodeGroup, findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.brainstem_terms import get_brainstem_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, Tapered, ConeBaseProgression, CylinderCentralPath
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.scaffoldpackage import ScaffoldPackage


class MeshType_3d_brainstem1(Scaffold_base):
    """
    Generates a tapered cylinder for the brainstem based on solid cylinder mesh, with variable numbers of elements in major, minor and length directions. Regions of the brainstem are annotated.
    """

    centralPathDefaultScaffoldPackages = {
        'Brainstem 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues( # dimensional.
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [ [0.0,0.0,-8.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-7.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-6.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-5.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-4.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-3.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-2.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0,-1.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ],
                    [ [0.0,0.0, 0.0], [0.0,0.0,1.0], [1.0,0.0,0.0], [0.0,0.0,0.0], [-0.0,1.0,0.0], [0.0,0.0,0.0] ]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_brainstem_term('medulla oblongata')[0],
                    'ontId': get_brainstem_term('medulla oblongata')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2',
                    'name': get_brainstem_term('pons')[0],
                    'ontId': get_brainstem_term('pons')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3',
                    'name': get_brainstem_term('midbrain')[0],
                    'ontId': get_brainstem_term('midbrain')[1]
                }]
        }),
        'Cat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues( # dimensional.
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [ [0.0,1.0,-5.0], [0.0,0.0,4.5], [5.0,0.0,0.0], [1.0,0.0,0.0], [0.0,2.4,0.0], [0.0,2.2,0.0] ],
                    [ [0.0,1.0,-0.5], [0.0,0.0,4.5], [6.0,0.0,0.0], [1.0,0.0,0.0], [0.0,4.0,0.0], [0.0,1.1,0.0] ],
                    [ [0.0,1.0, 4.0], [0.0,0.0,4.5], [7.0,0.0,0.0], [1.0,0.0,0.0], [0.0,4.5,0.0], [0.0,0.8,0.0] ],
                    [ [0.0,1.0, 8.5], [0.0,0.0,4.5], [8.0,0.0,0.0], [1.0,0.0,0.0], [0.0,5.5,0.0], [0.0,0.8,0.0] ],
                    [ [0.0,1.0,13.0], [0.0,0.0,4.5], [9.0,0.0,0.0], [1.0,0.0,0.0], [0.0,6.0,0.0], [0.0,0.2,0.0] ]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_brainstem_term('medulla oblongata')[0],
                    'ontId': get_brainstem_term('medulla oblongata')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2',
                    'name': get_brainstem_term('pons')[0],
                    'ontId': get_brainstem_term('pons')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3',
                    'name': get_brainstem_term('midbrain')[0],
                    'ontId': get_brainstem_term('midbrain')[1]
                }]
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 6
            },
            'meshEdits': exnodeStringFromNodeValues(  # dimensional.
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [0.0,18.9,-50.7], [0.0,-3.2,15.4], [ 7.8,-0.2,-0.0], [-0.1, 0.0,-0.1], [ 0.2, 6.1,1.3], [ 0.1, 0.5, 1.4] ],
                [ [0.0,14.6,-36.3], [0.0,-5.4,13.3], [ 8.5, 0.0,-0.0], [ 1.4, 0.3, 0.1], [ 0.0, 7.2,2.9], [-0.3, 1.8, 1.9] ],
                [ [0.0, 8.3,-24.1], [0.0,-6.4,12.5], [10.6, 0.4, 0.2], [ 4.8, 0.6, 0.2], [-0.5, 9.7,5.0], [-0.4, 1.6, 0.6] ],
                [ [0.0, 1.8,-11.4], [0.0,-5.3,13.3], [18.2, 1.1, 0.4], [ 4.8,-0.2,-0.1], [-0.7,10.4,4.2], [ 0.2, 3.0,-0.2] ],
                [ [0.0,-2.3,  2.4], [0.0,-3.9,13.2], [20.2, 0.0, 0.0], [-0.4, 0.9, 0.1], [-0.0,15.7,4.6], [-0.9, 2.3,-0.6] ],
                [ [0.0,-6.0, 15.0], [0.0,-2.6,12.9], [17.6, 2.7, 0.6], [-0.5, 1.4, 0.0], [-2.5,15.2,3.1], [-1.1,-0.6,-2.1] ],
                [ [0.0,-7.5, 28.0], [0.0,-0.4,13.1], [19.2, 2.9, 0.1], [ 3.7,-1.1,-1.0], [-2.2,14.5,0.4], [ 1.6,-0.7,-3.2] ]

                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_brainstem_term('medulla oblongata')[0],
                    'ontId': get_brainstem_term('medulla oblongata')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_brainstem_term('pons')[0],
                    'ontId': get_brainstem_term('pons')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-6',
                    'name': get_brainstem_term('midbrain')[0],
                    'ontId': get_brainstem_term('midbrain')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 6
            },
            'meshEdits': exnodeStringFromNodeValues(  # dimensional.
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [ [0.0, 2.9,-5.4], [0.0, 0.5,1.7], [1.9,0.0,0.0], [-0.2,0.2,0.1], [0.0,1.3,-0.3], [-0.1,0.3, 1.2] ],
                    [ [0.0, 2.7,-3.5], [0.0,-0.9,2.0], [1.8,0.1,0.0], [ 0.0,0.0,0.0], [-0.1,1.4,0.6], [-0.0,0.0, 0.7] ],
                    [ [0.0, 1.1,-1.8], [0.0,-1.4,1.7], [2.0,0.0,0.0], [ 0.5,0.0,0.0], [-0.0,1.3,1.1], [ 0.0,0.1, 0.2] ],
                    [ [0.0,-0.1,-0.2], [0.0,-1.1,1.5], [2.7,0.0,0.0], [ 0.3,0.0,0.0], [-0.0,1.5,1.1], [ 0.0,0.3,-0.1] ],
                    [ [0.0,-1.1, 1.2], [0.0,-0.7,1.6], [2.7,0.0,0.0], [-0.2,0.0,0.0], [-0.0,1.8,0.8], [ 0.0,0.1,-0.1] ],
                    [ [0.0,-1.5, 3.0], [0.0,-1.1,1.9], [2.2,0.0,0.0], [-0.1,0.0,0.0], [-0.0,1.7,1.0], [-0.0,-0.4,0.4] ],
                    [ [0.0,-3.3, 4.6], [0.0,-2.4,1.3], [2.7,0.0,0.1], [ 1.1,0.1,0.1], [-0.1,0.9,1.7], [-0.1,-1.2,0.9] ]
                    ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_brainstem_term('medulla oblongata')[0],
                    'ontId': get_brainstem_term('medulla oblongata')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_brainstem_term('pons')[0],
                    'ontId': get_brainstem_term('pons')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-6',
                    'name': get_brainstem_term('midbrain')[0],
                    'ontId': get_brainstem_term('midbrain')[1]
                }]
        })
    }

    @staticmethod
    def getName():
        return '3D Brainstem 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Brainstem 1',
            'Cat 1',
            'Human 1',
            'Rat 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if parameterSetName == 'Default':
            parameterSetName = 'Brainstem 1'

        if 'Brainstem 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Brainstem 1']

        elif 'Cat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Cat 1']

        elif 'Human 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']

        elif 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']

        centralPathOption1 = cls.centralPathDefaultScaffoldPackages['Brainstem 1']

        options = {
            'Base parameter set': parameterSetName,
            'Central path': copy.deepcopy(centralPathOption),
            'Brainstem path': copy.deepcopy(centralPathOption1),
            'Number of elements across major': 6,
            'Number of elements across minor': 6,
            'Number of elements along': 8,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major and minor': 1,
            'Refine number of elements along': 1
        }

        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            'Number of elements along',
            'Refine',
            'Refine number of elements across major and minor',
            'Refine number of elements along'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path' or optionName == 'Brainstem path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path' or optionName == 'Brainstem path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
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
        if optionName == 'Central path' or optionName == 'Brainstem path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False
        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1
        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1
        if options['Number of elements along'] < 2:
            options['Number of elements along'] = 2
        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isBrainstem = 'Brainstem 1' in parameterSetName
        isCat = 'Cat 1' in parameterSetName
        isHuman = 'Human 1' in parameterSetName
        isRat = 'Rat 1' in parameterSetName

        centralPath = options['Central path']
        brainstemPath = options['Brainstem path']
        elementsCountAcrossMajor = options['Number of elements across major']
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAlong = options['Number of elements along']

        # Cross section at Z axis
        halfBrainStem = False
        if halfBrainStem:
            elementsCountAcrossMajor //= 2
        elementsPerLayer = ((elementsCountAcrossMajor - 2) * elementsCountAcrossMinor) + (
                2 * (elementsCountAcrossMinor - 2))

        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        coordinates = findOrCreateFieldCoordinates(fm)
        mesh = fm.findMeshByDimension(3)

        # Annotation groups
        brainstemGroup = AnnotationGroup(region, get_brainstem_term('brainstem'))
        brainstemMeshGroup = brainstemGroup.getMeshGroup(mesh)
        midbrainGroup = AnnotationGroup(region, get_brainstem_term('midbrain'))
        midbrainMeshGroup = midbrainGroup.getMeshGroup(mesh)
        ponsGroup = AnnotationGroup(region, get_brainstem_term('pons'))
        ponsMeshGroup = ponsGroup.getMeshGroup(mesh)
        medullaGroup = AnnotationGroup(region, get_brainstem_term('medulla oblongata'))
        medullaMeshGroup = medullaGroup.getMeshGroup(mesh)
        annotationGroups = [brainstemGroup, midbrainGroup, ponsGroup, medullaGroup]

        annotationGroupAlong = [[brainstemGroup, midbrainGroup],
                                [brainstemGroup, ponsGroup],
                                [brainstemGroup, medullaGroup]]

        # point markers
        centralCanal = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('central canal of spinal cord'))
        cerebralAqueduct = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('cerebral aqueduct'))
        foramenCaecum = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('foramen caecum of medulla oblongata'))

        dorsalMidCaudalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('brainstem dorsal midline caudal point'))
        ventralMidCaudalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('brainstem ventral midline caudal point'))
        dorsalMidCranGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('brainstem dorsal midline cranial point'))
        ventralMidCranGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('brainstem ventral midline cranial point'))
        dorsalMidMedullaPonsJunction = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term ('brainstem dorsal midline pons-medulla junction'))
        ventralMidMedullaPonsJunction = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_brainstem_term('brainstem ventral midline pons-medulla junction'))
        dorsalMidMidbrainPonsJunction = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                      get_brainstem_term('brainstem dorsal midline midbrain-pons junction'))
        ventralMidMidbrainPonsJunction = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_brainstem_term('brainstem ventral midline midbrain-pons junction'))

        #######################
        # CREATE MAIN BODY MESH
        #######################
        brainTermsAlong = [None, 'medulla oblongata', 'pons', 'midbrain']
        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if not halfBrainStem else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        # Body coordinates
        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong, annotationGroup=brainTermsAlong)
        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor,
                            centre=[0.0, 0.0, 0.0],
                            alongAxis=cylinderCentralPath.alongAxis[0],
                            majorAxis=cylinderCentralPath.majorAxis[0],
                            minorRadius=cylinderCentralPath.minorRadii[0])

        cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath,
                                 useCrossDerivatives=False)

        brainstem_coordinates = findOrCreateFieldCoordinates(fm, name="brainstem coordinates")

        # Brain coordinates
        tmp_region = region.createRegion()
        tmp_fm = tmp_region.getFieldmodule()
        tmp_brainstem_coordinates = findOrCreateFieldCoordinates(tmp_fm, name="brainstem coordinates")

        cylinderCentralPath1 = CylinderCentralPath(tmp_region, brainstemPath, elementsCountAlong, annotationGroup=brainTermsAlong)

        base1 = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor,
                             centre=[0.0, 0.0, 0.0],
                             alongAxis=cylinderCentralPath1.alongAxis[0],
                             majorAxis=cylinderCentralPath1.majorAxis[0],
                             minorRadius=cylinderCentralPath1.minorRadii[0])

        cylinder2 = CylinderMesh(tmp_fm, tmp_brainstem_coordinates, elementsCountAlong, base1,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath1,
                                 useCrossDerivatives=False)

        # Write two coordinates
        sir = tmp_region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        tmp_region.write(sir)
        result, buffer = srm.getBuffer()

        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemoryBuffer(buffer)
        fielditer = fm.createFielditerator()
        field = fielditer.next()
        region.read(sir)

        while field.isValid():
            field = fielditer.next()
        del srm
        del sir
        del tmp_fm
        del tmp_brainstem_coordinates
        del tmp_region

        # Annotating groups
        iRegionBoundaries = [int(6 * elementsCountAlong / 15), int(13 * elementsCountAlong / 15)]
        for elementIdentifier in range(1, mesh.getSize() + 1):
            element = mesh.findElementByIdentifier(elementIdentifier)
            brainstemMeshGroup.addElement(element)
            if elementIdentifier > (iRegionBoundaries[-1] * elementsPerLayer):
                midbrainMeshGroup.addElement(element)
            elif (elementIdentifier > (iRegionBoundaries[0] * elementsPerLayer)) and (
                    elementIdentifier <= (iRegionBoundaries[-1] * elementsPerLayer)):
                ponsMeshGroup.addElement(element)
            else:
                medullaMeshGroup.addElement(element)

        ################
        # point markers
        ################
        pointMarkers = [
            {"group" : dorsalMidCaudalGroup, "marker_brainstem_coordinates" : [ 0.0, 1.0, -8.0 ]},
            {"group": ventralMidCaudalGroup, "marker_brainstem_coordinates": [0.0, -1.0, -8.0]},
            {"group": dorsalMidCranGroup, "marker_brainstem_coordinates": [0.0, 1.0, 0.0]},
            {"group": ventralMidCranGroup, "marker_brainstem_coordinates": [0.0, -1.0, 0.0]},
            {"group": dorsalMidMedullaPonsJunction, "marker_brainstem_coordinates": [0.0, 1.0, -5.0]},
            {"group": ventralMidMedullaPonsJunction, "marker_brainstem_coordinates": [0.0, -1.0, -5.0]},
            {"group": dorsalMidMidbrainPonsJunction, "marker_brainstem_coordinates": [0.0, 1.0, -2.0]},
            {"group": ventralMidMidbrainPonsJunction, "marker_brainstem_coordinates": [0.0, -1.0, -2.0]},
        ]

        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerBrainstemCoordinates = findOrCreateFieldCoordinates(fm, name="marker_body_coordinates")
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)
        markerTemplateInternal.defineField(markerBrainstemCoordinates)

        cache = fm.createFieldcache()

        brainstemNodesetGroup = brainstemGroup.getNodesetGroup(nodes)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        findMarkerLocation = fm.createFieldFindMeshLocation(markerBrainstemCoordinates, brainstem_coordinates, mesh)
        findMarkerLocation.setSearchMode(FieldFindMeshLocation.SEARCH_MODE_EXACT)

        for pointMarker in pointMarkers:
            group = pointMarker["group"]
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            cache.setNode(markerPoint)

            markerBrainstemCoordinates.assignReal(cache, pointMarker["marker_brainstem_coordinates"])
            markerName.assignString(cache, group.getName())

            element, xi = findMarkerLocation.evaluateMeshLocation(cache, 3)
            markerLocation.assignMeshLocation(cache, element, xi)
            group.getNodesetGroup(nodes).addNode(markerPoint)
            brainstemNodesetGroup.addNode(markerPoint)
            nodeIdentifier += 1

        return annotationGroups

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major and minor']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the 1D mesh.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for ventral-level elements.
        New point annotation groups are appended to this list.
        """
        # create  groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi1 = fm.createFieldOr(
            fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)),
            fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)))
        is_exterior_face_xi3 = fm.createFieldOr(fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)), fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)))

        # external regions
        namelist = ['brainstem', 'midbrain', 'medulla oblongata', 'pons']
        for subregion in namelist:
            subGroup = AnnotationGroup(region, get_brainstem_term(subregion))
            issub = subGroup.getFieldElementGroup(mesh2d)
            is_subface_ext = fm.createFieldOr(fm.createFieldAnd(issub, is_exterior_face_xi1), fm.createFieldAnd(issub, is_exterior_face_xi3))
            subFaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_brainstem_term(subregion + ' exterior'))
            subFaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_subface_ext)

def createCranialNerveEmergentMarkers(region, mesh, coordinatesName):
    # create marker points for locations the cranial nerves emerge from brainstem mesh, based on the USF cat brainstem data.
    # return element xi
    # use findMeshLocation to find the elementxi in an arbitrary mesh of given number of elements.

    if coordinatesName == "bodyCoordinates":
        # brainstem_coordinates: the left-side nerves
        nerveDict = {'OCULOMOTOR_left': [-0.13912257342955267, -0.5345161733750351, -0.7374762051676923],
                     'TROCHLEAR_left': [-0.13148279719950992, 0.4218745504359067, -0.7375838988856348],
                     'TRIGEMINAL_left': [-0.7605971693047597, -0.4025791045292648, -0.6862730212268676],
                     'ABDUCENS_left': [-0.19517975766630574, -0.6252563181242173, -0.8205128205130072],
                     'FACIAL_left': [-0.5824675040481234, -0.3554448371502354, -0.24509655553058302],
                     'VESTIBULOCOCHLEAR_left': [-0.6147505791411602, -0.32790803815838, -0.24509655403515848],
                     'GLOSSOPHARYNGEAL_left': [-0.7307312460087607, -0.2576952819028721, -0.39215539053073717],
                     'VAGUS_left': [-0.6741855912315219, -0.25981298010131126, -0.24509655277992023],
                     'ACCESSORY_cranialRoot_left': [-0.6741855912315219, -0.25981298010131126, -0.24509655277992023],
                     'HYPOGLOSSAL_left': [-0.044776303107883636, -0.5027870527016534, -0.10510117079651562]
                     }

        rightDict = {}
        for key in nerveDict.keys():
            nerveName = key.split('_')[0]+'_right'
            xyz = [-1*nerveDict[key][0], nerveDict[key][1], nerveDict[key][2]]
            rightDict.update({nerveName:xyz})
        nerveDict.update(rightDict)

        nerveNames = list(nerveDict.keys())

        # add to data coordinates
        markerNameField = 'marker_name'
        fm = region.getFieldmodule()
        cache = fm.createFieldcache()
        datapoints = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        data_coordinates = findOrCreateFieldCoordinates(fm, "data_coordinates")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        dnodetemplate = datapoints.createNodetemplate()
        dnodetemplate.defineField(data_coordinates)
        dnodetemplate.setValueNumberOfVersions(data_coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        dnodetemplate.defineField(markerName)
        dnodeIdentifier = 1
        for nerveName in nerveNames:
            node = datapoints.createNode(dnodeIdentifier, dnodetemplate)
            cache.setNode(node)
            addEnd = nerveDict[nerveName].copy()
            data_coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, addEnd)
            markerName.assignString(cache, nerveName)
            dnodeIdentifier += 1

        # find element-xi for these data_coordinates
        dataNamesField = fm.findFieldByName(markerNameField)
        coordinates = findOrCreateFieldCoordinates(fm, coordinatesName)
        found_mesh_location = fm.createFieldFindMeshLocation(data_coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)
        xi_projected_data = {}
        nodeIter = datapoints.createNodeiterator()
        node = nodeIter.next()
        while node.isValid():
            cache.setNode(node)
            element, xi = found_mesh_location.evaluateMeshLocation(cache, 3)
            marker_name = dataNamesField.evaluateString(cache)
            if element.isValid():
                addProjection = {marker_name:{"elementID": element.getIdentifier(), "xi": xi,"nodeID": node.getIdentifier()}}
                xi_projected_data.update(addProjection)
            node = nodeIter.next()

        result = datapoints.destroyAllNodes()

    else:
        xi_projected_data = {}

    return xi_projected_data

