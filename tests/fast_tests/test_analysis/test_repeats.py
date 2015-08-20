'''
Tests Repeats analysis class.

'''

from nose.tools import assert_equal
from pymbt import analysis, DNA


def test_find_repeats():
    input_sequence = DNA('atgatgccccgatagtagtagtag')
    expected = [('ATG', 2), ('GTA', 3), ('GAT', 2), ('AGT', 3), ('CCC', 2),
                ('TAG', 4)]

    output = analysis.repeats(input_sequence, 3)
    assert_equal(output, expected)
