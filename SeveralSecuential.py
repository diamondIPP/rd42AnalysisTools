#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from TestAreas import TestAreas

class SeveralSecuential:
	def __init__(self):
		self.runs = [25202, 25203, 25205, 25206, 25207, 25208, 25209, 25210]
		self.out_dir = '/home/sandiego/data/2018/output/no_sat_50_6'
		self.config_path = '/home/sandiego/rd42AnalysisTools/AreaSettings/2018_10/config_50um_test1.ini'

		for run in self.runs:
			testArea = TestAreas(self.config_path, run)
			testArea.DoNoiseStudiesDifferentBuffers('good')

if __name__ == '__main__':
	t = SeveralSecuential()
