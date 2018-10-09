
import pytest

from SICER2.src.reads_to_bins import add_reads_to_dict, files_to_bins


args = {"bin_size": 200, "fragment_size": 150, "drop_duplicates": True}


def test_files_to_bins():

    d = files_to_bins(["tests/test.bed"], args)
    print(d)

    assert list(d["chr7", "+"]) == [20246600, 91135000]


def test_reads_to_dict():

    d = add_reads_to_dict("tests/test.bed")
    print(d)

    assert list(d["chr7", "+"]) == [20246668, 20246668, 91135110]
