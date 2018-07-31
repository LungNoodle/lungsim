import os
import unittest

from aether.diagnostics import set_diagnostics_on
from aether.geometry import define_node_geometry_2d
from aether.arrays import check_node_xyz_2d

# Look to see if the 'TEST_RESOURCES_DIR' is set otherwise fallback to a path
# relative to this file.
if 'TEST_RESOURCES_DIR' in os.environ:
    resources_dir = os.environ['TEST_RESOURCES_DIR']
else:
    here = os.path.abspath(os.path.dirname(__file__))
    resources_dir = os.path.join(here, 'resources')


class GeometryTestCase(unittest.TestCase):

    def test_read_square(self):
        set_diagnostics_on(False)
        define_node_geometry_2d(os.path.join(resources_dir, 'square.ipnode'))
        value = check_node_xyz_2d(1, 1, 10)
        self.assertEqual(10, value)


if __name__ == '__main__':
    unittest.main()
