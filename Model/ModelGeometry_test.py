import unittest
import ModelGeometry
import numpy as np
class TestGeometry(unittest.TestCase):
    def test_nodedistance(self):
        dz = ModelGeometry.nodedistance(np.array([0,1]), 2)
        self.assertEqual(dz, 1)

    def test_choosegeometry(self):
        [nz, Z, coord] = ModelGeometry.choose_geometry('FieldScale0.5m')
        self.assertEqual(nz,101)
        self.assertEqual(Z, 0.5)
        self.assertEqual(sum(coord), sum(np.linspace(0,Z,nz)))

if __name__ == '__main__':
    unittest.main()