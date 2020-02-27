# Third party
import pytest

# First party
from propy.AAIndex import GetAAIndex1, GetAAIndex23


@pytest.mark.skip(reason="Currently fails on travis-ci.org")
def test_main():
    # init(path='.')
    # grep('volume')
    # x = get('KRIW790103')
    # print(x)
    # print(x.get('W'))
    temp1 = GetAAIndex1("KRIW790103")
    print(len(temp1))

    temp2 = GetAAIndex23("TANS760101")
    print(len(temp2))
    temp2 = GetAAIndex23("GRAR740104")
    print(len(temp2))
