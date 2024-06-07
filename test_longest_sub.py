import pytest
from longest_sub import longest_contiguous, all_longest, run_large_test_parallel
from Bio.Seq import Seq


@pytest.mark.parametrize("seq, stop_codon, gap, case_sensitive, expected, expected_exception", [
    (Seq("a"), "", "", True, {"a": 1}, None),
    (Seq("aabbbbcc"), "", "", True, {"b": 4}, None),
    (Seq("aabb"), "", "", True, {"a": 2, "b": 2}, None),
    (Seq(""), "", "", True, {}, None),
    (Seq("a***b"), "*", "", True, {"*": 3}, None),
    (Seq("ab--c"), "", "-", True, {"-": 2}, None),
    (Seq("abcAAA"), "", "", True, {"A": 3}, None),
    (Seq("aaaaBBbbBBA"), "", "", True, {"a": 4}, None), # case sensitive 
    (Seq("aaaaBBbbBBA"), "", "", False, {"b": 6}, None), # not case sensitive 


    # cases expected to raise ValueError
    (Seq("AB-456$"), "", "", True, None, ValueError), # invalid characters 
    (Seq("*78 **"), "", "", True, None, ValueError), # space in sequence
    ("a", "", "", True, None, ValueError), # none Seq type
    ("", "", "", True, None, ValueError), # empty argument 
    (Seq("AUGAAA"), "AUG", "", True, None, ValueError), # stop codon longer than 1
    (Seq("AUGAAA"), "", "GAP", True, None, ValueError), # gap longer than 1
])


def test_longest_contiguous(seq, stop_codon, gap, case_sensitive, expected, expected_exception):
    if expected_exception:
        with pytest.raises(expected_exception):
            longest_contiguous(seq, stop_codon=stop_codon, gap=gap, case_sensitive=case_sensitive)
    else:
        output = longest_contiguous(seq, stop_codon=stop_codon, gap=gap, case_sensitive=case_sensitive)
        assert output == expected



@pytest.mark.parametrize("records, expected", [
    ([{"a": 4}], {"a": 4}),
    ([{"a": 4}, {"b": 6, "a": 2}, {"b": 10}], {"a": 4, "b": 10}),
    ([], {}),
    ([{}], {}),
])
def test_all_longest(records, expected):
    output = all_longest(records)

    assert output == expected

def run_large_test_parallel(records, expected):
    output = run_large_test_parallel(records)

    assert output == expected


def test_integration():
    output = all_longest(longest_contiguous(seq) for seq in [Seq("aaa"), Seq("aabbbbcc"), Seq("aabb")])

    assert output == {"a": 3, "b": 4}