"""
Generates a solid ellipsoid of hexahedral elements.
"""
import math
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.ellipsoidmesh import EllipsoidMesh, EllipsoidSurfaceD3Mode
from scaffoldmaker.utils.meshgeneratedata import MeshGenerateData
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_ellipsoid1(Scaffold_base):
    """
    Generates a solid ellipsoid of hexahedral elements.
    """

    @classmethod
    def getName(cls):
        return "3D Ellipsoid 1"

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            "Numbers of elements across axes": [4, 6, 8],
            "Numbers of shell, transition elements": [0, 1],
            "Axes lengths": [1.0, 1.5, 2.0],
            "Axes shell thicknesses": [0.2, 0.2, 0.2],
            "Axis 2 x-rotation degrees": 0.0,
            "Axis 3 x-rotation degrees": 90.0,
            "Advanced n-way derivative factor": 0.6,
            "Advanced surface D3 mode": EllipsoidSurfaceD3Mode.SURFACE_NORMAL.value,
            "Core": True,
            "Core shell scaling mode": 1,
            "Use linear through shell": False,
            "Refine": False,
            "Refine number of elements": 4,
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Numbers of elements across axes",
            "Numbers of shell, transition elements",
            "Axes lengths",
            "Axes shell thicknesses",
            "Axis 2 x-rotation degrees",
            "Axis 3 x-rotation degrees",
            "Core",
            "Core shell scaling mode",
            "Advanced n-way derivative factor",
            "Advanced surface D3 mode",
            # "Use linear through shell",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependent_changes = False

        max_rim_count = None
        axes_numbers = options["Numbers of elements across axes"]
        count = len(axes_numbers)
        if count < 3:
            for i in range(3 - count):
                axes_numbers.append(axes_numbers[-1])
        elif count > 3:
            del axes_numbers[3:]
        for i, number in enumerate(axes_numbers):
            if number < 4:
                axes_numbers[i] = 4
            elif number % 2:
                axes_numbers[i] += 1
            transition_count = (axes_numbers[i] // 2) - 1
            if (max_rim_count is None) or (transition_count < max_rim_count):
                max_rim_count = transition_count

        shell_transition_counts = options["Numbers of shell, transition elements"]
        count = len(axes_numbers)
        if count < 2:
            shell_transition_counts[1] = 1
        elif count > 2:
            del shell_transition_counts[2:]
        if shell_transition_counts[0] > max_rim_count - 1:
            shell_transition_counts[0] = max_rim_count - 1
            dependent_changes = True
        if shell_transition_counts[1] < 1:
            shell_transition_counts[1] = 1
        elif (shell_transition_counts[1] + shell_transition_counts[0]) > max_rim_count:
            shell_transition_counts[1] = max_rim_count - shell_transition_counts[0]
            dependent_changes = True

        axes_lengths = options["Axes lengths"]
        count = len(axes_lengths)
        if count < 3:
            for i in range(3 - count):
                axes_lengths.append(axes_lengths[-1])
        elif count > 3:
            del axes_lengths[3:]
        for i, length in enumerate(axes_lengths):
            if length <= 0.0:
                axes_lengths[i] = 1.0

        axes_thicknesses = options["Axes shell thicknesses"]
        count = len(axes_thicknesses)
        if count < 3:
            for i in range(3 - count):
                axes_thicknesses.append(axes_thicknesses[-1])
        elif count > 3:
            del axes_thicknesses[3:]
        for i, thickness in enumerate(axes_thicknesses):
            if thickness <= 0.0:
                axes_thicknesses[i] = axes_lengths[i] * 0.2

        if options["Core shell scaling mode"] not in (1, 2):
            options["Core shell scaling mode"] = 1
        if options["Advanced n-way derivative factor"] < 0.1:
            options["Advanced n-way derivative factor"] = 0.1
        elif options["Advanced n-way derivative factor"] > 1.0:
            options["Advanced n-way derivative factor"] = 1.0

        try:
            mode = EllipsoidSurfaceD3Mode(options["Advanced surface D3 mode"])
        except ValueError:
            options["Advanced surface D3 mode"] = EllipsoidSurfaceD3Mode.SURFACE_NORMAL.value

        for key in [
            "Refine number of elements"
        ]:
            if options[key] < 1:
                options[key] = 1

        return dependent_changes

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: empty list of AnnotationGroup, None
        """
        element_counts = options["Numbers of elements across axes"]
        shell_element_count, transition_element_count = options["Numbers of shell, transition elements"]
        axes_lengths = options["Axes lengths"]
        axes_shell_thicknesses = options["Axes shell thicknesses"]

        axis2_x_rotation_radians = math.radians(options["Axis 2 x-rotation degrees"])
        axis3_x_rotation_radians = math.radians(options["Axis 3 x-rotation degrees"])
        core = options["Core"]
        core_shell_scaling_mode = options["Core shell scaling mode"]
        nway_d_factor = options["Advanced n-way derivative factor"]
        surface_d3_mode = EllipsoidSurfaceD3Mode(options["Advanced surface D3 mode"])

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)

        ellipsoid = EllipsoidMesh(element_counts, shell_element_count, transition_element_count, core,
                                  core_shell_scaling_mode=core_shell_scaling_mode)

        left_group = AnnotationGroup(region, ("left", ""))
        right_group = AnnotationGroup(region, ("right", ""))
        back_group = AnnotationGroup(region, ("back", ""))
        front_group = AnnotationGroup(region, ("front", ""))
        bottom_group = AnnotationGroup(region, ("bottom", ""))
        top_group = AnnotationGroup(region, ("top", ""))
        annotation_groups = [left_group, right_group, back_group, front_group, bottom_group, top_group]
        octant_group_lists = []
        for octant in range(8):
            octant_group_list = []
            octant_group_list.append((right_group if (octant & 1) else left_group).getGroup())
            octant_group_list.append((front_group if (octant & 2) else back_group).getGroup())
            octant_group_list.append((top_group if (octant & 4) else bottom_group).getGroup())
            octant_group_lists.append(octant_group_list)
        ellipsoid.set_octant_group_lists(octant_group_lists)

        if core:
            box_group = AnnotationGroup(region, ("box", ""))
            transition_group = AnnotationGroup(region, ("transition", ""))
            annotation_groups += [box_group, transition_group]
            ellipsoid.set_box_transition_groups(box_group.getGroup(), transition_group.getGroup())

        ellipsoid.build(axes_lengths, axis2_x_rotation_radians, axis3_x_rotation_radians, axes_shell_thicknesses,
                        nway_d_factor=nway_d_factor, surface_d3_mode=surface_d3_mode)
        generate_data = MeshGenerateData(region, meshDimension=(2 if ((shell_element_count == 0) and not core) else 3))
        ellipsoid.generate_mesh(generate_data)

        return annotation_groups, None

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCount = options["Refine number of elements"]
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)
