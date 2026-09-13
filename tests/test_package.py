import re

import openabel


def test_version_is_semver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", openabel.__version__)


def test_public_api_is_abel_only():
    assert openabel.__all__ == ["Abel"]
