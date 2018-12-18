

import openAbel as oa
from nose.tools import *


@raises(NotImplementedError)
def test_methodNotImplemented():
    oa.Abel(10, -1, 0., 1., method = -1)


@raises(NotImplementedError)
def test_methodNotImplemented():
    oa.Abel(10, -1, 0., 1., method = 5)
