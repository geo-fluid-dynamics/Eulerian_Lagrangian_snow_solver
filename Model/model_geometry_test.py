import unittest
import model_geometry
import numpy as np
class TestGeometry(unittest.TestCase):
    def test_nodedistance(self):
        dz = model_geometry.node_distance(np.array([0,1]), 2)
        self.assertEqual(dz, 1)

    def test_choosegeometry(self):
        [nz,Z, coord] = model_geometry.choose_geometry('FieldScale0.5m')
        self.assertEqual(nz,101)
        self.assertEqual(Z, 0.5)
        self.assertEqual(sum(coord), sum(np.linspace(0,Z,nz)))

if __name__ == '__main__':
    unittest.main()