import os
from nose.tools import assert_equal, assert_raises, assert_true, assert_false
import coral as cr


def test_construction():
    plasmid_path = os.path.join(os.path.dirname(__file__), 'gibson_test.fasta')
    f1_path = os.path.join(os.path.dirname(__file__), 'fragment_1.fasta')
    f2_path = os.path.join(os.path.dirname(__file__), 'fragment_2.fasta')
    f3_path = os.path.join(os.path.dirname(__file__), 'fragment_3.fasta')
    f3_linear_path = os.path.join(os.path.dirname(__file__),
                                  'fragment_3_linear.fasta')
    plasmid = cr.io.read_dna(plasmid_path).circularize()
    f1 = cr.io.read_dna(f1_path)
    f2 = cr.io.read_dna(f2_path)
    f3 = cr.io.read_dna(f3_path)
    f3_linear = cr.io.read_dna(f3_linear_path)

    gibsoned_circular = cr.reaction.gibson([f1, f2, f3])
    gibsoned_linear = cr.reaction.gibson([f1, f2, f3_linear], linear=True)
#
    expected_length = len(plasmid)
    gibsoned_circular_length = len(gibsoned_circular)
    gibsoned_linear_length = len(gibsoned_linear)
    assert_equal(gibsoned_circular_length, expected_length)
    assert_equal(gibsoned_linear_length, expected_length)
    assert_true(gibsoned_circular.circular)
    assert_false(gibsoned_linear.circular)
    assert(plasmid.is_rotation(gibsoned_circular))
    try:
        assert_equal(str(plasmid), str(gibsoned_linear))
    except AssertionError:
        assert_equal(str(plasmid), str(gibsoned_linear.flip()))

    # Should fail with circular input
    assert_raises(ValueError, cr.reaction.gibson, [f1.circularize()])
    # Should fail if compatible end can't be found
    assert_raises(Exception, cr.reaction.gibson, [f1, f3[50:-50]], linear=True)
    normal = [f1, f2, f3]
    rotated = [f1, f2, f3.reverse_complement()]
    # Gibson should work regardless of fragment orientation
    gibsoned = cr.reaction.gibson(normal)
    assert_true(gibsoned.is_rotation(cr.reaction.gibson(rotated)))
    # A redundant fragment shouldn't affect the outcome
    assert_equal(cr.reaction.gibson([f1, f2, f3]),
                 cr.reaction.gibson([f1, f2, f2, f3]))
    # A fragment that can't circularize should raise a ValueError
    assert_raises(ValueError, cr.reaction.gibson, [f1, f2, f3[:-80]])
    # But should still work fine as a linear fragment


def test_annotations():
    pass
