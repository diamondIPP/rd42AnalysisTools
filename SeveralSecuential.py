#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from optparse import OptionParser
from TestAreas import TestAreas
from Utils import *

class SeveralSecuential:
	def __init__(self):
		self.runs = [25202, 25203, 25205, 25206, 25207, 25208, 25209, 25210]
		self.out_dir = '/home/sandiego/data/2018/output/no_sat_50_6'
		self.config_path = '/home/sandiego/rd42AnalysisTools/AreaSettings/2018_10/config_50um_test1.ini'
		self.ta = {run: TestAreas(self.config_path, run) for run in self.runs}


	def ConfigureRuns(self):
		for run, testArea in self.ta.iteritems():
			testArea.SetTransparentGrid()
			testArea.SetTest()
			if testArea.trans_grid.loaded_pickle:
				testArea.trans_grid.LoadPickle()
			else:
				testArea.trans_grid.FindPickleValues()
				ExitMessage('Run it again to load the generated pickle :)', os.EX_OK)

	def DoRuns(self):
		for run, testArea in self.ta.iteritems():
			testArea.DoNoiseStudiesDifferentBuffers('good')

if __name__ == '__main__':
	t = SeveralSecuential()
	t.ConfigureRuns()
	t.DoRuns()