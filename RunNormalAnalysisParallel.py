#!/usr/bin/env python
import time
import subprocess as subp
import multiprocessing as mp
import os, shutil
import ipdb
from optparse import OptionParser
from copy import deepcopy
import numpy as np
from collections import OrderedDict
from ParallelManager import ParallelManager


class RunNormalAnalysisParallel:
	def __init__(self, runlist, num_cores=2, force=False, verb=True):
		self.verb = verb
		self.runlist = runlist
		self.force = force
		self.num_cores = num_cores if num_cores <= int(mp.cpu_count()/2.0) or force else int(mp.cpu_count()/2.0)
		self.num_runs = 0
		self.settings_list = []
		self.job_chunks = []
		self.analysis_processes = {}
		self.workind_dir = os.getcwd()
		self.queue = {}
		self.queue_running = {}
		self.queue_showing = {}
		self.queue_runs = {}
		self.runs_dic_completed = {}
		self.runs_dic_running = {}
		if not os.path.isfile(self.runlist):
			print 'File', self.runlist, 'does not exist. Exiting'
			exit(os.EX_CONFIG)
		print 'Starting parallel analysis using runlist', self.runlist, 'and using', self.num_cores, 'cores simultaneously'
		self.ReadRunList()
		self.parallelManager = ParallelManager()
		self.RunParallelAnalysis()
		print 'Runs completed. Exiting'

	def ReadRunList(self):
		with open(self.runlist, 'r') as rl:
			lines = rl.readlines()
			self.settings_list = [line.replace('\n', '') for line in lines if os.path.isfile(line.replace('\n', ''))]
			self.num_runs = len(self.settings_list)
			self.job_chunks = [self.settings_list[i:i + self.num_cores] for i in xrange(0, self.num_runs, self.num_cores)]
			self.num_cores = min(self.num_cores, self.num_runs)
			# print 'Jobs are grouped as following:', self.job_chunks

	def RunParallelAnalysis(self):
		working_dir = os.getcwd()
		options = [['-w', working_dir, '-s', os.path.abspath(jobi), '--normal'] for jobi in self.settings_list]
		if not self.verb:
			options = [option + ['-q'] for option in options]

		self.parallelManager.SetVariables(working_dir=working_dir, runlist=self.settings_list, exec_command='rd42Analysis.py', options=options, num_cores=self.num_cores, force=self.force, verb=self.verb)
		self.parallelManager.RunParallelAnalysis()

def main():
	parser = OptionParser()
	parser.add_option('-r', '--runlist', dest='runlist', type='string', help='File containing a list of the RunSettings for each run')
	parser.add_option('-n', '--numcores', dest='numcores', type='int', default=2, help='number of runs to execute in parallel')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true', help='force to use the specified number of cores')
	parser.add_option('-q', '--quiet', dest='quiet', default=False, action='store_true', help='enables quiet mode: no verbose')

	(options, args) = parser.parse_args()
	runlist = options.runlist
	num = options.numcores
	force = options.force
	verb = not bool(options.quiet)

	pp = RunNormalAnalysisParallel(runlist=runlist, num_cores=num, force=force, verb=verb)


if __name__ == '__main__':
	main()
