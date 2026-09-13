import re

import openAbel


def test_version_isSemver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", openAbel.__version__)


def test_publicApi_isAbelOnly():
    assert openAbel.__all__ == ["Abel"]
