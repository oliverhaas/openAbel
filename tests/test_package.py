import re

import openabel


def test_version_isSemver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", openabel.__version__)


def test_publicApi_isAbelOnly():
    assert openabel.__all__ == ["Abel"]
