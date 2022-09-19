import bedmaker
import pytest
import os


class TestSmoke:
    """
    Testing input parser
    """

    def test_correct(self):
        assert True


cor_dir = "data/bed/correct"
bed_files_correct = [
    os.path.join(cor_dir, d) for d in os.listdir(cor_dir)
]
incor_dir = "data/bed/incorrect"
bed_files_incorrect = [
    os.path.join(incor_dir, d) for d in os.listdir(incor_dir)
]

class TestBedqc:
    """
    Testing input parser
    """

    @pytest.mark.parametrize("bedfile", bed_files_correct)
    def test_correct(self, bedfile, tmpdir):
        error_list = bedmaker.bedqc.run_bedqc(bedfile, tmpdir)
        assert len(error_list) == 0

    @pytest.mark.parametrize("bedfile", bed_files_incorrect)
    def test_incorrect(self, bedfile, tmpdir):
        error_list = bedmaker.bedqc.run_bedqc(bedfile, tmpdir)
        assert len(error_list) > 0
