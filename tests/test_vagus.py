from cmlibs.maths.vectorops import mult
from cmlibs.utils.zinc.field import (
    find_or_create_field_coordinates, find_or_create_field_finite_element, find_or_create_field_group,
    find_or_create_field_stored_string)
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.group import mesh_group_to_identifier_ranges, nodeset_group_to_identifier_ranges
from cmlibs.zinc.context import Context
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK

from scaffoldmaker.annotation.annotationgroup import findAnnotationGroupByName
from scaffoldmaker.annotation.vagus_terms import vagus_branch_terms, vagus_marker_terms
from scaffoldmaker.meshtypes.meshtype_3d_nerve1 import MeshType_3d_nerve1, get_left_vagus_marker_locations_list
from scaffoldmaker.utils.interpolation import get_curve_from_points, getCubicHermiteCurvesLength
from scaffoldmaker.utils.read_vagus_data import VagusInputData
from testutils import assertAlmostEqualList, check_annotation_term_ids

import math
import os
import unittest


here = os.path.abspath(os.path.dirname(__file__))


def reorder_vagus_test_data1(testcase, region):
    """
    Break up and reorder vagus test data to test trunk ordering code.
    """
    fieldmodule = region.getFieldmodule()
    mesh1d = fieldmodule.findMeshByDimension(1)
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    trunk_group = fieldmodule.findFieldByName("left vagus X nerve trunk").castGroup()
    coordinates = fieldmodule.findFieldByName("coordinates")
    IDENTIFIER_OFFSET = 1000
    UNUSED_IDENTIFIER = 100000
    with ChangeManager(fieldmodule):
        for element_identifier in range(1, 101):
            element = mesh1d.findElementByIdentifier(element_identifier)
            eft = element.getElementfieldtemplate(coordinates, -1)
            local_node_count = eft.getNumberOfLocalNodes()
            node_identifiers = [element.getNode(eft, ln).getIdentifier() for ln in range(1, local_node_count + 1)]
            node_identifiers.reverse()
            testcase.assertEqual(RESULT_OK, element.setNodesByIdentifier(eft, node_identifiers))
        # don't renumber nodes 46, 60 as they're the start of a branch
        # currently rely on them being the first identifier in that branch for a workaround which includes it
        for node_identifier in range(1, 40):
            other_node_identifier = 101 - node_identifier
            node = nodes.findNodeByIdentifier(node_identifier)
            other_node = nodes.findNodeByIdentifier(other_node_identifier)
            testcase.assertEqual(RESULT_OK, other_node.setIdentifier(UNUSED_IDENTIFIER))
            testcase.assertEqual(RESULT_OK, node.setIdentifier(other_node_identifier))
            testcase.assertEqual(RESULT_OK, other_node.setIdentifier(node_identifier))
        for element_identifier in range(76, 101):
            other_element_identifier = 151 - element_identifier
            element = mesh1d.findElementByIdentifier(element_identifier)
            other_element = mesh1d.findElementByIdentifier(other_element_identifier)
            testcase.assertEqual(RESULT_OK, other_element.setIdentifier(UNUSED_IDENTIFIER))
            testcase.assertEqual(RESULT_OK, element.setIdentifier(other_element_identifier))
            testcase.assertEqual(RESULT_OK, other_element.setIdentifier(element_identifier))
        # for node_identifier in range(126, 151):
        #     node = nodes.findNodeByIdentifier(node_identifier)
        #     testcase.assertEqual(RESULT_OK, node.setIdentifier(node_identifier + IDENTIFIER_OFFSET))
        for element_identifier in range(101, 104):
            mesh1d.destroyElement(mesh1d.findElementByIdentifier(element_identifier))


class VagusScaffoldTestCase(unittest.TestCase):


    def test_vagus_terms(self):
        """
        Test nomenclature of the vagus terms. 
        """
        for term_ids in (vagus_branch_terms + vagus_marker_terms):
            self.assertTrue(check_annotation_term_ids(term_ids), "Invalid primary term id or order not UBERON < ILX < FMA for vagus annotation term ids " + str(term_ids)) 

    def test_input_vagus_data(self):
        """
        Reading vagus input data
        """

        data_file = os.path.join(here, "resources", "vagus_test_data1.exf")

        # test with original ordered, and reordered vagus data file
        for i in range(2):
            context = Context("Test")
            base_region = context.getDefaultRegion()
            region = base_region.createChild('vagus')
            assert(region.isValid())
            data_region = region.getParent().createChild('data')
            assert(data_region.isValid())
            result = data_region.readFile(data_file)
            assert result == RESULT_OK
            if i == 1:
                reorder_vagus_test_data1(self, data_region)

            data_fieldmodule = data_region.getFieldmodule()
            data_mesh1d = data_fieldmodule.findMeshByDimension(1)
            data_nodes = data_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            trunk_group = data_fieldmodule.findFieldByName("left vagus X nerve trunk").castGroup()
            data_coordinates = data_fieldmodule.findFieldByName("coordinates")
            data_trunk_mesh_group = trunk_group.getMeshGroup(data_mesh1d)
            data_trunk_nodeset_group = trunk_group.getNodesetGroup(data_nodes)
            mesh_ranges = mesh_group_to_identifier_ranges(data_trunk_mesh_group)
            nodeset_ranges = nodeset_group_to_identifier_ranges(data_trunk_nodeset_group)
            if i == 0:
                self.assertEqual([[1, 200]], mesh_ranges)
                self.assertEqual([[1, 201]], nodeset_ranges)
                self.assertEqual(200, data_trunk_mesh_group.getSize())
                expected_element_info = {
                    1: [1, 2],
                    2: [2, 3],
                    51: [51, 52],
                    101: [101, 102],
                    102: [102, 103],
                    103: [103, 104]
                }
            else:
                self.assertEqual([[1, 100], [104, 200]], mesh_ranges)
                # self.assertEqual([[1, 125], [151, 201], [1126, 1150]], nodeset_ranges)
                self.assertEqual(197, data_trunk_mesh_group.getSize())
                expected_element_info = {
                    1: [99, 100],
                    2: [98, 99],
                    51: [101, 1],
                    101: None,
                    102: None,
                    103: None
                }
            for element_id, expected_node_ids in expected_element_info.items():
                element = data_mesh1d.findElementByIdentifier(element_id)
                eft = element.getElementfieldtemplate(data_coordinates, -1)
                if expected_node_ids:
                    node_ids = [element.getNode(eft, n + 1).getIdentifier() for n in range(eft.getNumberOfLocalNodes())]
                    self.assertEqual(expected_node_ids, node_ids)
                else:
                    self.assertFalse(eft.isValid())

            vagus_data = VagusInputData(data_region)
            self.assertEqual(vagus_data.get_side_label(), 'left')

            marker_data = vagus_data.get_level_markers()
            self.assertEqual(len(marker_data), 4)
            self.assertTrue('left level of superior border of the clavicle on the vagus nerve' in marker_data)

            orientation_data = vagus_data.get_orientation_data()
            self.assertEqual(len(orientation_data), 8)
            expected_orientation_info = [
                ("orientation left", 1),
                ("orientation left anterior", 2),
                ("orientation anterior", 12),
                ("orientation right anterior", 1),
                ("orientation right", 1),
                ("orientation right posterior", 1),
                ("orientation posterior", 1),
                ("orientation left posterior", 1),
            ]
            for orientation_direction_name, expected_count in expected_orientation_info:
                self.assertEqual(len(orientation_data[orientation_direction_name]), expected_count)

            trunk_group_name = vagus_data.get_trunk_group_name()
            self.assertEqual(trunk_group_name, "left vagus nerve")
            trunk_coordinates = vagus_data.get_trunk_coordinates()
            self.assertEqual(len(trunk_coordinates), 201)
            annotation_term_map = vagus_data.get_annotation_term_map()
            self.assertTrue(trunk_group_name in annotation_term_map)
            # self.assertEqual(annotation_term_map[trunk_group_name], 'http://purl.obolibrary.org/obo/UBERON_0035020')

            # do a simple fit to the trunk data coordinates to check trunk ordering is working
            trunk_data_coordinates = vagus_data.get_trunk_coordinates()
            px = [e[0] for e in trunk_data_coordinates]
            self.assertEqual(201, len(px))
            bx, bd1 = get_curve_from_points(px, number_of_elements=10)
            length = getCubicHermiteCurvesLength(bx, bd1)
            self.assertAlmostEqual(31726.825262197974, length, delta=1.0E-3)

            branch_coordinates_data = vagus_data.get_branch_coordinates_data()
            branch_sequences_data = vagus_data.get_branch_sequences_data()
            self.assertEqual(len(branch_coordinates_data), 4)
            self.assertTrue("left superior laryngeal nerve" in branch_coordinates_data)
            self.assertEqual(len(branch_coordinates_data["left superior laryngeal nerve"]), 44)
            self.assertEqual(len(branch_sequences_data["left superior laryngeal nerve"]), 2)
            self.assertTrue("left A branch of superior laryngeal nerve" in branch_coordinates_data)
            self.assertEqual(len(branch_coordinates_data["left A branch of superior laryngeal nerve"]), 22)
            left_thoracic_cardiopulmonary_branches = (
                "left A thoracic cardiopulmonary branch of vagus nerve",
                "left B thoracic cardiopulmonary branch of vagus nerve")
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertTrue(branch_name in branch_coordinates_data)

            branch_parents = vagus_data.get_branch_parent_map()
            self.assertEqual(branch_parents["left superior laryngeal nerve"], trunk_group_name)
            self.assertEqual(branch_parents["left A branch of superior laryngeal nerve"], "left superior laryngeal nerve")
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertEqual(branch_parents[branch_name], trunk_group_name)

            branch_common_groups = vagus_data.get_branch_common_group_map()
            self.assertEqual(len(branch_common_groups), 2)
            self.assertTrue("left branch of superior laryngeal nerve" in branch_common_groups)
            self.assertTrue("left A branch of superior laryngeal nerve" in \
                            branch_common_groups["left branch of superior laryngeal nerve"])
            self.assertTrue("left thoracic cardiopulmonary branch of vagus nerve" in branch_common_groups)
            for branch_name in left_thoracic_cardiopulmonary_branches:
                self.assertTrue(branch_name in branch_common_groups["left thoracic cardiopulmonary branch of vagus nerve"])

    def test_no_input_file(self):
        """
        No input file.
        """

        scaffold = MeshType_3d_nerve1
        options = scaffold.getDefaultOptions("Human Right Vagus 1")

        context = Context("Test")
        base_region = context.getDefaultRegion()
        region = base_region.createChild('vagus')

        # check that it doesn't crash with  no input file
        annotation_groups = scaffold.generateBaseMesh(region, options)[0]
        self.assertEqual(annotation_groups, [])

    def test_vagus_nerve_1(self):
        """
        Test creation of vagus nerve scaffold with simple, synthetic data similar to REVA data.
        """
        scaffold = MeshType_3d_nerve1
        scaffoldname = MeshType_3d_nerve1.getName()
        self.assertEqual(scaffoldname, '3D Nerve 1')
        parameterSetNames = scaffold.getParameterSetNames()
        self.assertEqual(parameterSetNames, ['Default', 'Human Left Vagus 1', 'Human Right Vagus 1'])
        options = scaffold.getDefaultOptions("Human Left Vagus 1")
        self.assertEqual(len(options), 8)
        self.assertEqual(options.get('Base parameter set'), 'Human Left Vagus 1')
        self.assertEqual(options.get('Number of elements along the trunk pre-fit'), 30)
        self.assertEqual(options.get('Number of elements along the trunk'), 50)
        self.assertEqual(options.get('Trunk proportion'), 1.0)
        self.assertEqual(options.get('Trunk fit number of iterations'), 5)
        self.assertEqual(options.get('Default anterior direction'), [0.0, 1.0, 0.0])
        self.assertEqual(options.get('Default trunk diameter'), 3.0)
        self.assertEqual(options.get('Branch diameter trunk proportion'), 0.5)
        # change options to make test fast and consistent, with minor effect on result:
        options['Number of elements along the trunk pre-fit'] = 10
        options['Number of elements along the trunk'] = 25
        options['Trunk fit number of iterations'] = 2
        options['Default trunk diameter'] = 300.0

        # test with original ordered, and reordered vagus data file
        for i in range(2):

            context = Context("Test")
            root_region = context.getDefaultRegion()
            region = root_region.createChild('vagus')
            data_region = root_region.createChild('data')
            data_file = os.path.join(here, "resources", "vagus_test_data1.exf")
            self.assertEqual(data_region.readFile(data_file), RESULT_OK)
            if i == 1:
                reorder_vagus_test_data1(self, data_region)

            # create segment groups dividing the data approximately in thirds over the x-span of the trunk
            data_fieldmodule = data_region.getFieldmodule()
            data_nodes = data_fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            data_coordinates = data_fieldmodule.findFieldByName("coordinates")
            with ChangeManager(data_fieldmodule):
                data_x = data_fieldmodule.createFieldComponent(data_coordinates, 1)
                conditions = [
                    data_fieldmodule.createFieldLessThan(data_x, data_fieldmodule.createFieldConstant(10000.0)),
                    None,
                    data_fieldmodule.createFieldGreaterThan(data_x, data_fieldmodule.createFieldConstant(20000.0))
                    ]
                conditions[1] = data_fieldmodule.createFieldNot(
                    data_fieldmodule.createFieldOr(conditions[0], conditions[2]))
                for s in range(3):
                    segment_group = data_fieldmodule.createFieldGroup()
                    segment_group.setName("segment" + str(s + 1) + ".exf")
                    segment_group.setManaged(True)
                    segment_nodeset_group = segment_group.createNodesetGroup(data_nodes)
                    segment_nodeset_group.addNodesConditional(conditions[s])
                del conditions
                del data_x

            # check annotation groups
            annotation_groups, nerve_metadata = scaffold.generateMesh(region, options)
            self.assertEqual(len(annotation_groups), 20)
            metadata = nerve_metadata.getMetadata()["vagus nerve"]
            TOL = 1.0E-6
            expected_metadata = {
                'segments': {
                    'segment1.exf': {'minimum vagus coordinate': 0.060890363163301214,
                                     'maximum vagus coordinate': 0.245177232602641},
                    'segment2.exf': {'minimum vagus coordinate': 0.24778486671128517,
                                     'maximum vagus coordinate': 0.40913081735379395},
                    'segment3.exf': {'minimum vagus coordinate': 0.41170171261995313,
                                     'maximum vagus coordinate': 0.5737279929406741}
                },
                'trunk centroid fit error rms': 2.9700127717877277,
                'trunk centroid fit error max': 11.107680110311745,
                'trunk radius fit error rms': 1.3533423007380196,
                'trunk radius fit error max': 5.428778876520852,
                'trunk twist angle fit error degrees rms': 3.9145909417227405,
                'trunk twist angle fit error degrees max': 9.787350215289262}
            self.assertEqual(len(metadata), len(expected_metadata))
            for key, value in metadata.items():
                expected_value = expected_metadata[key]
                if key == 'segments':
                    self.assertEqual(len(value), len(expected_value))
                    for segment_name, vagus_coordinate_range in value.items():
                        expected_vagus_coordinate_range = expected_value[segment_name]
                        for range_key, range_value in vagus_coordinate_range.items():
                            self.assertAlmostEqual(range_value, expected_vagus_coordinate_range[range_key], delta=TOL)
                else:
                    self.assertAlmostEqual(value, expected_value, delta=TOL)

            # (term_id, parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3,
            #  expected_surface_area, expected_volume)
            expected_group_info = {
                'left vagus nerve': (
                    'http://uri.interlex.org/base/ilx_0785628', None, 25,
                    [-1207.6880100059116, -6475.308801043239, -58.4223048224834],
                    [2279.4644089925127, -983.6137056743332, 109.78753490210283],
                    [42.63617203080673, 274.1720386906218, 1571.1413726765452],
                    251553677.33283976,
                    33870558978.020557),
                'left superior laryngeal nerve': (
                    'http://uri.interlex.org/base/ilx_0788780', 'left vagus nerve', 4,
                    [5923.104657597035, -4450.2479197707235, -196.91175665569313],
                    [-1473.6650516759184, 858.0807042974038, 37.61890734384133],
                    [29.505599722138413, 24.943359756955033, 586.8845479660235],
                    14603873.1104719,
                    895690578.4083588),
                'left A branch of superior laryngeal nerve': (
                    'http://uri.interlex.org/base/ilx_0795823', 'left superior laryngeal nerve', 2,
                    [5105.456364262517, -1456.2684055690102, 0.18793093373062675],
                    [-1289.5812951072821, 381.4601337342459, 17.49393056176476],
                    [2.9906263178490917, -3.6461869711379222, 299.9629335014554],
                    4696615.984152875,
                    236649715.5919552),
                'left A thoracic cardiopulmonary branch of vagus nerve': (
                    'http://uri.interlex.org/base/ilx_0794192', 'left vagus nerve', 2,
                    [20637.123232118378, -2947.094130818922, -608.0143068866589],
                    [99.3811573594038, -1713.8817535655432, -61.05879554434705],
                    [-8.82246083735481, 11.986772318349722, -350.8208201352616],
                    6221084.957871486,
                    330616841.6518934),
                'left B thoracic cardiopulmonary branch of vagus nerve': (
                    'http://uri.interlex.org/base/ilx_0794193', 'left vagus nerve', 1,
                    [22164.372546340623, -3219.413785808221, -620.4335804281038],
                    [1775.1656728388946, 1620.626138221386, -217.2367728432068],
                    [2.313286113960203, 43.732207832213135, 345.15346649012787],
                    4676303.367303241,
                    269749902.83920854)
            }
            groups_count = len(expected_group_info)

            # check all meshes and groups are created and of appropriate size
            fieldmodule = region.getFieldmodule()
            fieldcache = fieldmodule.createFieldcache()
            coordinates = fieldmodule.findFieldByName("coordinates").castFiniteElement()
            self.assertTrue(coordinates.isValid())
            self.assertEqual(RESULT_OK, fieldmodule.defineAllFaces())
            mesh3d = fieldmodule.findMeshByDimension(3)
            expected_elements_count = 34
            self.assertEqual(expected_elements_count, mesh3d.getSize())
            mesh2d = fieldmodule.findMeshByDimension(2)
            # groups_count + 1 due to one group having 2 branches
            self.assertEqual(expected_elements_count * 9 + groups_count + 1, mesh2d.getSize())
            mesh1d = fieldmodule.findMeshByDimension(1)
            self.assertEqual(expected_elements_count * 17 + (groups_count + 1) * 8, mesh1d.getSize())
            nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
            self.assertEqual(expected_elements_count + 1 + 8, nodes.getSize())  # including 6 marker points
            datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
            self.assertEqual(0, datapoints.getSize())
            meshDerivative1 = mesh3d.getChartDifferentialoperator(order=1, term=1)
            meshDerivative3 = mesh3d.getChartDifferentialoperator(order=1, term=3)
            XTOL = 0.01  # coordinates and derivatives
            LTOL = 0.001  # length
            STOL = 0.1  # surface area
            VTOL = 1.0  # volume
            MTOL = 1.0E-7  # material coordinate
            trunk_group_name = "left vagus trunk"
            one = fieldmodule.createFieldConstant(1.0)
            for group_name in expected_group_info.keys():
                term_id, parent_group_name, expected_elements_count, expected_start_x, expected_start_d1, expected_start_d3, \
                    expected_surface_area, expected_volume = expected_group_info[group_name]
                group = fieldmodule.findFieldByName(group_name).castGroup()
                annotation_group = findAnnotationGroupByName(annotation_groups, group_name)
                self.assertEqual(term_id, annotation_group.getId())
                mesh_group3d = group.getMeshGroup(mesh3d)
                elements_count = mesh_group3d.getSize()
                expected_face_count = expected_elements_count * 9 + 1
                expected_line_count = expected_elements_count * 17 + 8
                expected_node_count = expected_elements_count + (2 if parent_group_name else 1)
                if group_name == 'left superior laryngeal nerve':
                    # there are 2 branches with this name
                    expected_face_count += 1
                    expected_line_count += 8
                    expected_node_count += 1
                self.assertEqual(expected_elements_count, elements_count)
                mesh_group2d = group.getMeshGroup(mesh2d)
                self.assertEqual(expected_face_count, mesh_group2d.getSize(), msg=group_name)
                mesh_group1d = group.getMeshGroup(mesh1d)
                self.assertEqual(expected_line_count, mesh_group1d.getSize(), msg=group_name)
                nodeset_group = group.getNodesetGroup(nodes)
                self.assertEqual(expected_node_count, nodeset_group.getSize(), msg=group_name)
                branch_of_branch = False
                if parent_group_name:
                    # check first 2 nodes are in parent nodeset group
                    parent_group = fieldmodule.findFieldByName(parent_group_name).castGroup()
                    parent_nodeset_group = parent_group.getNodesetGroup(nodes)
                    nodeiterator = nodeset_group.createNodeiterator()
                    for n in range(2):
                        node = nodeiterator.next()
                        self.assertTrue(parent_nodeset_group.containsNode(node))
                    branch_of_branch = parent_group_name != trunk_group_name
                element = mesh_group3d.createElementiterator().next()
                self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5]))
                result, start_x = coordinates.evaluateReal(fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d1 = coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                TOL = 10.0 * XTOL if branch_of_branch else XTOL
                assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
                assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
                assertAlmostEqualList(self, start_d3, expected_start_d3, delta=TOL)
                # note surface area is merely sum of all surfaces including epineurium, box and interior surfaces
                surface_area_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group2d)
                surface_area_field.setNumbersOfPoints(4)
                volume_field = fieldmodule.createFieldMeshIntegral(one, coordinates, mesh_group3d)
                volume_field.setNumbersOfPoints(3)
                fieldcache.clearLocation()  # otherwise integrates over element only
                result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                result, volume = volume_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                self.assertAlmostEqual(expected_surface_area, surface_area, delta=350.0 if branch_of_branch else STOL)
                self.assertAlmostEqual(expected_volume, volume, delta=20000.0 if branch_of_branch else VTOL)

            # check sampled trunk d3 for orientation and radius fit, all at element centre
            xi_centre = [0.5, 0.5, 0.5]
            # (element_identifier, expected_d3)
            expected_d3_info = [
                (2, [-28.137582921730427, 259.5994524928721, 1174.0350511301763]),
                (4, [-514.5998529487936, 672.3468661596938, 819.2553817278631]),
                (6, [-29.808187017002297, 642.3244458693601, 197.31879904198632]),
                (8, [-23.562152113000536, 202.73422030358648, 666.4595633643032]),
                (10, [-29.248031392808684, -281.10525266406444, 642.7009578428102]),
                (12, [-243.16037712446342, -571.2413263097393, 242.1937017613861]),
                (14, [-0.36150950965947004, -475.284647464952, 117.36033062680346]),
                (16, [-2.511317112796206, -469.3593933416148, 106.48467255534705])]
            for element_identifier, expected_d3 in expected_d3_info:
                element = mesh3d.findElementByIdentifier(element_identifier)
                self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, xi_centre))
                result, d3 = coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                assertAlmostEqualList(self, d3, expected_d3, delta=XTOL)

            # check volume of trunk, surface area of epineurium, length of centroids, coordinates and straight coordinates
            straight_coordinates = fieldmodule.findFieldByName("straight coordinates").castFiniteElement()
            self.assertTrue(straight_coordinates.isValid())
            STOL = 100.0
            LTOL = 0.5
            for coordinate_field in (coordinates, straight_coordinates):
                group = fieldmodule.findFieldByName("left vagus nerve").castGroup()
                mesh_group3d = group.getMeshGroup(mesh3d)
                volume_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group3d)
                volume_field.setNumbersOfPoints(3)
                fieldcache.clearLocation()
                result, volume = volume_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                expected_volume = 33870558978.020557 if (coordinate_field is coordinates) else 33873441076.58173
                self.assertAlmostEqual(expected_volume, volume, delta=STOL)
                expected_elements_count = 34
                group = fieldmodule.findFieldByName("epineurium").castGroup()
                mesh_group2d = group.getMeshGroup(mesh2d)
                self.assertEqual(expected_elements_count * 4, mesh_group2d.getSize())
                surface_area_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group2d)
                surface_area_field.setNumbersOfPoints(4)
                fieldcache.clearLocation()
                result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                expected_surface_area = 74241915.15945683 if (coordinate_field is coordinates) else 74385368.0658539
                self.assertAlmostEqual(expected_surface_area, surface_area, delta=STOL)
                group = fieldmodule.findFieldByName("vagus centroid").castGroup()
                mesh_group1d = group.getMeshGroup(mesh1d)
                self.assertEqual(expected_elements_count, mesh_group1d.getSize())
                length_field = fieldmodule.createFieldMeshIntegral(one, coordinate_field, mesh_group1d)
                length_field.setNumbersOfPoints(4)
                result, length = length_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                self.assertAlmostEqual(77872.46497432186, length, delta=LTOL)

            # check all markers are added
            marker_group = fieldmodule.findFieldByName("marker").castGroup()
            marker_nodes = marker_group.getNodesetGroup(nodes)
            self.assertEqual(8, marker_nodes.getSize())
            marker_name_field = fieldmodule.findFieldByName("marker_name").castStoredString()
            marker_location_field = fieldmodule.findFieldByName("marker_location").castStoredMeshLocation()
            expected_marker_info = get_left_vagus_marker_locations_list()
            expected_elements_count = 25
            expected_marker_node_identifier = 10001
            node_iter = marker_nodes.createNodeiterator()
            node = node_iter.next()
            for expected_marker_name, expected_material_coordinate3 in expected_marker_info.items():
                self.assertEqual(expected_marker_node_identifier, node.getIdentifier())
                fieldcache.setNode(node)
                marker_name = marker_name_field.evaluateString(fieldcache)
                self.assertEqual(expected_marker_name, marker_name)
                element, xi = marker_location_field.evaluateMeshLocation(fieldcache, 3)
                material_coordinate3 = (element.getIdentifier() - 1 + xi[0]) / expected_elements_count
                self.assertAlmostEqual(expected_material_coordinate3, material_coordinate3, delta=MTOL)
                expected_marker_node_identifier += 1
                node = node_iter.next()

            # check vagus material coordinates
            vagus_coordinates = fieldmodule.findFieldByName("vagus coordinates").castFiniteElement()
            self.assertTrue(vagus_coordinates.isValid())
            # (expected_start_x, expected_start_d1, expected_start_d3, expected_surface_area, expected_volume)
            expected_group_material_info = {
                'left vagus nerve': (
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.04],
                    [0.0, 0.012, 0.0],
                    0.07044881379783888,
                    0.00014399999999999916),
                'left superior laryngeal nerve': (
                    [0.00048701207606156654, 0.00016977995042870763, 0.1332411665654758],
                    [0.012580581637176687, 0.011551948251027854, -0.00593608421881044],
                    [-0.004039218464322834, 0.004386752607670381, -8.224378638307939e-05],
                    0.0029560000433909083,
                    2.957815257008607e-06),
                'left A branch of superior laryngeal nerve': (
                    [0.028289902111450684, 0.024975964871614482, 0.11922843754888073],
                    [-0.009456996982071633, -0.01317644600853053, -0.017368559424983603],
                    [-0.0042481671245299946, 0.0041499700558325844, -0.0008953060937929858],
                    0.0015939616370540632,
                    1.5192676822737654e-06),
                'left A thoracic cardiopulmonary branch of vagus nerve': (
                    [-0.00023180370608450157, -8.872684747093166e-06, 0.38043566433142256],
                    [-0.026494164248939335, 0.01090954289884798, 0.0061236730593394115],
                    [-0.0022767025683330944, -0.005549330555353137, -3.464648527395919e-05],
                    0.0020623089007186924,
                    2.1093381596184995e-06),
                'left B thoracic cardiopulmonary branch of vagus nerve': (
                    [0.0005186352499216176, -0.0009142039153784858, 0.40544076796886575],
                    [0.0233995801517628, -0.026607181414355724, 0.02017770849343402],
                    [0.004505535415037674, 0.00396216809721049, 1.8417935306325184e-05],
                    0.0014331096986807375,
                    1.470024527726565e-06)}
            XTOL = 1.0E-4  # coordinates and derivatives
            STOL = 1.0E-5  # surface area
            VTOL = 1.0E-8  # volume
            for group_name in expected_group_info.keys():
                expected_start_x, expected_start_d1, expected_start_d3, expected_surface_area, expected_volume = \
                    expected_group_material_info[group_name]
                group = fieldmodule.findFieldByName(group_name).castGroup()
                mesh_group3d = group.getMeshGroup(mesh3d)
                mesh_group2d = group.getMeshGroup(mesh2d)
                element = mesh_group3d.createElementiterator().next()
                self.assertEqual(RESULT_OK, fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5]))
                result, start_x = vagus_coordinates.evaluateReal(fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d1 = vagus_coordinates.evaluateDerivative(meshDerivative1, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                result, start_d3 = vagus_coordinates.evaluateDerivative(meshDerivative3, fieldcache, 3)
                self.assertEqual(RESULT_OK, result)
                branch_of_branch = (group_name == 'left A branch of superior laryngeal nerve')
                TOL = (10.0 * XTOL) if branch_of_branch else XTOL
                assertAlmostEqualList(self, start_x, expected_start_x, delta=TOL)
                assertAlmostEqualList(self, start_d1, expected_start_d1, delta=TOL)
                assertAlmostEqualList(self, start_d3, expected_start_d3, delta=TOL)
                # note surface area is merely sum of all surfaces including epineurium, box and interior surfaces
                surface_area_field = fieldmodule.createFieldMeshIntegral(one, vagus_coordinates, mesh_group2d)
                surface_area_field.setNumbersOfPoints(4)
                volume_field = fieldmodule.createFieldMeshIntegral(one, vagus_coordinates, mesh_group3d)
                volume_field.setNumbersOfPoints(3)
                fieldcache.clearLocation()  # otherwise integrates over element only
                result, surface_area = surface_area_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                result, volume = volume_field.evaluateReal(fieldcache, 1)
                self.assertEqual(result, RESULT_OK)
                self.assertAlmostEqual(expected_surface_area, surface_area, delta=2.0E-6 if branch_of_branch else STOL)
                self.assertAlmostEqual(expected_volume, volume, delta=2.0E-9 if branch_of_branch else VTOL)

            # test combined groups
            branch_common_map = {
                'left branch of superior laryngeal nerve': [
                    'left A branch of superior laryngeal nerve'],
                'left thoracic cardiopulmonary branch of vagus nerve': [
                    'left A thoracic cardiopulmonary branch of vagus nerve',
                    'left B thoracic cardiopulmonary branch of vagus nerve']}
            for common_branch_name, variant_branch_name_list in branch_common_map.items():
                common_mesh_group = fieldmodule.findFieldByName(common_branch_name).castGroup().getMeshGroup(mesh3d)
                sum_variant_group_sizes = 0
                for variant_branch_name in variant_branch_name_list:
                    variant_mesh_group = fieldmodule.findFieldByName(variant_branch_name).castGroup().getMeshGroup(mesh3d)
                    sum_variant_group_sizes += variant_mesh_group.getSize()
                    elem_iter = variant_mesh_group.createElementiterator()
                    element = elem_iter.next()
                    while element.isValid():
                        self.assertTrue(common_mesh_group.containsElement(element))
                        element = elem_iter.next()
                self.assertEqual(common_mesh_group.getSize(), sum_variant_group_sizes)

            # test sizes of cervical and thoracic parts
            expected_section_info = {
                "left cervical vagus nerve": ("http://uri.interlex.org/base/ilx_0794142", 8),
                "left thoracic vagus nerve": ("http://uri.interlex.org/base/ilx_0787543", 17)}
            for section_name, info in expected_section_info.items():
                expected_id, expected_mesh_size = info
                annotation_group = findAnnotationGroupByName(annotation_groups, section_name)
                self.assertEqual(expected_id, annotation_group.getId())
                self.assertEqual(expected_mesh_size, annotation_group.getMeshGroup(mesh3d).getSize())

    def test_arc_vagus(self):
        """
        Test creation of a vagus nerve scaffold following a half circle so longitudinal curvature
        effects are seen.
        """
        context = Context("Test")
        root_region = context.getDefaultRegion()
        region = root_region.createChild('vagus')
        data_region = root_region.createChild('data')
        generate_arc_vagus_data(data_region)

        scaffold = MeshType_3d_nerve1
        options = scaffold.getDefaultOptions('Human Left Vagus 1')
        elements_count = 16
        options['Number of elements along the trunk pre-fit'] = elements_count
        options['Number of elements along the trunk'] = elements_count
        options['Trunk fit number of iterations'] = 2

        annotation_groups, nerve_metadata = scaffold.generateMesh(region, options)
        self.assertEqual(14, len(annotation_groups))
        fit_metadata = nerve_metadata.getMetadata()['vagus nerve']
        self.assertAlmostEqual(fit_metadata['trunk centroid fit error rms'], 0.0, delta=1.0E-4)
        self.assertAlmostEqual(fit_metadata['trunk radius fit error rms'], 4.550429168566819e-05, delta=1.0E-12)
        self.assertAlmostEqual(fit_metadata['trunk twist angle fit error degrees rms'], 0.0, delta=0.002)
        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        mesh1d = fieldmodule.findMeshByDimension(1)
        coordinates = fieldmodule.findFieldByName('coordinates').castFiniteElement()
        self.assertTrue(coordinates.isValid())
        one = fieldmodule.createFieldConstant(1.0)
        centroid_group = fieldmodule.findFieldByName('vagus centroid').castGroup()
        self.assertTrue(centroid_group.isValid())
        centroid_mesh_group = centroid_group.getMeshGroup(mesh1d)
        self.assertEqual(elements_count, centroid_mesh_group.getSize())
        length_field = fieldmodule.createFieldMeshIntegral(one, coordinates, centroid_mesh_group)
        length_field.setNumbersOfPoints(4)
        result, length = length_field.evaluateReal(fieldcache, 1)
        self.assertEqual(result, RESULT_OK)
        self.assertAlmostEqual(math.pi, length, delta=1.0E-2)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        node = nodes.findNodeByIdentifier((elements_count // 2) + 1)
        fieldcache.setNode(node)
        XTOL = 1.0E-3
        expected_parameters = [
            [1.000000940622472, -6.338102288830355e-06, 0.0],
            [1.8235912738924405e-07, 0.19637019676125334, 0.0],
            [-0.07071023719107788, 6.566504166275671e-08, 0.07071111904345162],
            [6.81154153478958e-05, -0.01390580597932258, 6.812747932889646e-05],
            [0.07071111904342112, -6.566586059468453e-08, 0.07071023719110837],
            [6.814039311262761e-05, 0.013905979277005844, -6.812832897063537e-05]]
        i = 0
        for value_label in [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]:
            result, parameters = coordinates.getNodeParameters(fieldcache, -1, value_label, 1, 3)
            self.assertEqual(RESULT_OK, result)
            assertAlmostEqualList(self, expected_parameters[i], parameters, delta=XTOL)
            i += 1


def generate_arc_vagus_data(region):
    """
    Generate a semicircle of vagus data for testing centroid curvature terms in vagus nerve scaffold.
    The data has no branches.
    :param region: Region to create vagus data in.
    """
    fieldmodule = region.getFieldmodule()
    with (ChangeManager(fieldmodule)):
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh1d = fieldmodule.findMeshByDimension(1)
        coordinates = find_or_create_field_coordinates(fieldmodule, managed=True)
        radius = find_or_create_field_finite_element(fieldmodule, "radius", 1, managed=True)

        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.defineField(radius)
        elementtemplate = mesh1d.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        linear_basis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        eft = mesh1d.createElementfieldtemplate(linear_basis)
        elementtemplate.defineField(coordinates, -1, eft)
        elementtemplate.defineField(radius, -1, eft)

        trunk_group = find_or_create_field_group(fieldmodule, "left vagus X nerve trunk", managed=True)
        trunk_nodes = trunk_group.createNodesetGroup(nodes)
        trunk_mesh1d = trunk_group.createMeshGroup(mesh1d)

        half_elements_count = 100
        elements_count = 2 * half_elements_count
        fieldcache = fieldmodule.createFieldcache()
        r = 0.05
        for n in range(elements_count + 1):
            node = trunk_nodes.createNode(n + 1, nodetemplate)
            fieldcache.setNode(node)
            theta = 0.5 * math.pi * (n - half_elements_count) / half_elements_count
            x = [math.cos(theta), math.sin(theta), 0.0]
            coordinates.assignReal(fieldcache, x)
            radius.assignReal(fieldcache, r)
            if n > 0:
                element = trunk_mesh1d.createElement(n, elementtemplate)
                element.setNodesByIdentifier(eft, [n, n + 1])

        anterior_group = find_or_create_field_group(fieldmodule, "orientation anterior", managed=True)
        anterior_nodes = anterior_group.createNodesetGroup(nodes)
        node_identifier = elements_count + 2
        anterior_points_count = 20
        anterior_angle = math.pi / 4.0
        anterior_r = 1.0 + math.cos(math.pi / 4.0) * 1.5 * r
        anterior_z = math.sin(anterior_angle) * 1.5 * r
        for n in range(anterior_points_count):
            node = anterior_nodes.createNode(node_identifier, nodetemplate)
            fieldcache.setNode(node)
            theta = math.pi * (-0.5 + (n + 0.5) / anterior_points_count)
            x = [anterior_r * math.cos(theta), anterior_r * math.sin(theta), anterior_z]
            coordinates.assignReal(fieldcache, x)
            radius.assignReal(fieldcache, 0.5 * r)
            node_identifier += 1

        # add level marker data

        datapoints = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        marker_group = find_or_create_field_group(fieldmodule, "marker", managed=True)
        marker_datapoints = marker_group.createNodesetGroup(datapoints)
        marker_name = find_or_create_field_stored_string(fieldmodule, "marker_name", managed=True)
        nodetemplate_marker = datapoints.createNodetemplate()
        nodetemplate_marker.defineField(coordinates)
        nodetemplate_marker.defineField(marker_name)

        data_identifier = 1
        marker_info = get_left_vagus_marker_locations_list()
        for name, material_coordinate3 in marker_info.items():
            node = marker_datapoints.createNode(data_identifier, nodetemplate_marker)
            fieldcache.setNode(node)
            theta = math.pi * (-0.5 + material_coordinate3)
            x = [math.cos(theta), math.sin(theta), 0.0]
            coordinates.assignReal(fieldcache, x)
            marker_name.assignString(fieldcache, name)
            data_identifier += 1


if __name__ == "__main__":
    unittest.main()
