import math
import os
import tempfile

from aether.exports import export_data_geometry
from aether.geometry import define_data_geometry


def read_exported_points(filename):
    points = []
    with open(filename, encoding="utf-8") as stream:
        lines = iter(stream)
        for line in lines:
            if line.lstrip().startswith("Node:"):
                points.append(tuple(map(float, next(lines).split())))
    return points


def check_points(actual, expected):
    assert len(actual) == len(expected)
    for actual_point, expected_point in zip(actual, expected):
        assert all(
            math.isclose(actual_value, expected_value, rel_tol=0.0, abs_tol=1.0e-10)
            for actual_value, expected_value in zip(actual_point, expected_point)
        )


expected = [(1.25, -2.5, 3.75), (-4.0, 5.5, -6.25)]

with tempfile.TemporaryDirectory() as directory:
    exdata_base = os.path.join(directory, "points")
    with open(exdata_base + ".exdata", "w", encoding="utf-8") as stream:
        stream.write(
            " Group name: test\n"
            " #Fields=2\n"
            " 1) marker, field, rectangular cartesian, #Components=1\n"
            "   value. Value index= 1, #Derivatives=0\n"
            " 2) coordinates, coordinate, rectangular cartesian, #Components=3\n"
            "   x. Value index= 2, #Derivatives=0\n"
            "   y. Value index= 3, #Derivatives=0\n"
            "   z. Value index= 4, #Derivatives=0\n"
            " Node: 101\n"
            "  9.0 1.25 -2.5 3.75\n"
            " Node: 205\n"
            "  8.0 -4.0 5.5 -6.25\n"
        )

    define_data_geometry(exdata_base)
    output_base = os.path.join(directory, "from-exdata")
    export_data_geometry(output_base, "test", 0)
    check_points(read_exported_points(output_base + ".exdata"), expected)

    ipdata_file = os.path.join(directory, "points.ipdata")
    with open(ipdata_file, "w", encoding="utf-8") as stream:
        stream.write(
            " converted from exdata\n"
            " 1 1.25 -2.5 3.75 1.0 1.0 1.0\n"
            " 2 -4.0 5.5 -6.25 1.0 1.0 1.0\n"
        )

    define_data_geometry(ipdata_file)
    output_base = os.path.join(directory, "from-ipdata")
    export_data_geometry(output_base, "test", 0)
    check_points(read_exported_points(output_base + ".exdata"), expected)

print("EXDATA and IPDATA point import succeeded.")
