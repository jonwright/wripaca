
import wripaca
import unittest

class test_det3( unittest.TestCase ):
    def test1(self):
        v = wripaca.new_doubleArray(9)
        for i,x in enumerate([1,0,0,0,1,0,0,0,1]):
            wripaca.doubleArray_setitem( v, i, x )
        assert wripaca.determinant3( v ) == 1.0


def testsuite():
    testSuite = unittest.TestSuite()
    testSuite.addTest( test_det3( "test1"))
    return testSuite

if __name__=="__main__":
    mysuite = testsuite()
    runner = unittest.TextTestRunner()
    runner.run( mysuite )   
