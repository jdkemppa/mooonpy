# -*- coding: utf-8 -*-

import pytest
from mooonpy import Thermospace, Path

class TestThermospace:

    def test_read(self):
        """Test logfile read"""
        log = Thermospace.basic_read(Path('../examples/thermospace/dummy_log.txt'),True)
        assert log.shape() == (30,7)

    def test_join(self):
        """Test joining two logfiles"""
        log = Thermospace.basic_read(Path('../examples/thermospace/dummy_log.txt'),True)
        restart = Thermospace.basic_read(Path('../examples/thermospace/dummy_restart_log.txt'))
        log.join_restart(restart)
        sect = log.sect(None,'first')
        assert sect.shape() == (31,7)

# if __name__ == '__main__':
#     pytest.main()
    # log = Thermospace.basic_read(Path('../examples/thermospace/dummy_log.txt'),True)
    # restart  = Thermospace.basic_read(Path('../examples/thermospace/dummy_restart_log.txt'))
    # log.join_restart(restart)
    # sect = log.sect(None,'first')