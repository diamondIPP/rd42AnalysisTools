#!/usr/bin/env python
import time
import subprocess as subp
import multiprocessing as mp
import os, shutil
import ipdb
from optparse import OptionParser
from copy import deepcopy
import numpy as np

class RunNormalAnalysisParallel:
	def __init__(self, runlist, num_cores=2, force=False):
		self.runlist = runlist
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
		self.all_completed = False
		if not os.path.isfile(self.runlist):
			print 'File', self.runlist, 'does not exist. Exiting'
			exit(os.EX_CONFIG)
		print 'Starting parallel analysis using runlist', self.runlist, 'and using', self.num_cores, 'cores simultaneously'
		self.ReadRunList()
		self.RunParallelAnalysis2()

	def ReadRunList(self):
		with open(self.runlist, 'r') as rl:
			lines = rl.readlines()
			self.settings_list = [line.split('\n')[0] for line in lines if os.path.isfile(line.split('\n')[0])]
			self.num_runs = len(self.settings_list)
			self.job_chunks = [self.settings_list[i:i + self.num_cores] for i in xrange(0, self.num_runs, self.num_cores)]
			self.num_cores = min(self.num_cores, self.num_runs)
			print 'Jobs are grouped as following:', self.job_chunks

	def RunParallelAnalysis(self):
		with open(os.devnull, 'w') as FNULL:
			for jobs in self.job_chunks:
				self.analysis_processes = []
				for it, run in enumerate(jobs):
					if it == len(jobs) - 1:
						print 'Showing output for run', run
						self.analysis_processes.append(subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(run), '--normal'], bufsize=-1, stdin=subp.PIPE, close_fds=True))
					else:
						self.analysis_processes.append(subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(run), '--normal', '-q'], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True))
				for job_i in xrange(len(self.analysis_processes)):
					while self.analysis_processes[job_i].poll() is None:
						time.sleep(5)
					self.CloseSubprocess(self.analysis_processes[job_i], stdin=True, stdout=False)
				print 'Done with', jobs

	def RunParallelAnalysis2(self):
		self.runs_dic_completed = {run: False for run in self.settings_list}
		self.runs_dic_running = {run: False for run in self.settings_list}
		self.queue = {c: None for c in xrange(self.num_cores)}
		self.queue_running = {c: False for c in xrange(self.num_cores)}
		self.queue_showing = {c: False for c in xrange(self.num_cores)}
		self.queue_runs = {c: None for c in xrange(self.num_cores)}
		first_time = True
		with open(os.devnull, 'w') as FNULL:
			while not np.array(self.runs_dic_completed.values(), '?').all():
				pos_run = np.bitwise_xor(self.runs_dic_completed.values(), self.runs_dic_running.values()).argmin()
				with self.settings_list[pos_run] as jobi:
					do_add_queue = not np.array(self.queue_running.values(), '?').all()
					if do_add_queue:
						pos_q, nfree = np.array(self.queue_running.values(), '?').argmin(), np.bitwise_not(self.queue_running.values()).sum()
						if nfree == 1 and not np.array(self.queue_showing.values(), '?').any():
							print 'Showing output for run', jobi
							self.queue[pos_q] = subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(jobi), '--normal'], bufsize=-1, stdin=subp.PIPE, close_fds=True)
							self.queue_showing[pos_q] = True
						else:
							self.queue[pos_q] = subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(jobi), '--normal'], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True)
						self.queue_running[pos_q] = True
						self.runs_dic_running[jobi] = True
						self.queue_runs[pos_q] = jobi
					if not first_time:
						temp = deepcopy(self.queue_running)
						for p, queue_p in temp.itervalues():
							if queue_p:
								if self.queue[p]:
									temp2 = self.queue[p]
									if temp2.poll():
										self.CloseSubprocess(self.queue[p], stdin=True, stdout=False)
										self.queue[p] = None
										self.queue_running[p] = False
										if self.queue_showing[p]:
											self.queue_showing[p] = False
										with self.queue_runs[p] as jobj:
											self.runs_dic_running[jobj] = False
											self.runs_dic_completed[jobj] = True
						time.sleep(3)
					else:
						first_time = not np.array(self.queue_running.values(), '?').all()

	def GetAvailablePos(self):
		pos = -1
		num_free = 0
		for p, elem in self.queue.iteritems():
			if not elem:
				if pos == -1:
					pos = p
				num_free += 1
		return pos, num_free

	def CloseSubprocess(self, p, stdin=False, stdout=False):
		pid = p.pid
		if stdin:
			p.stdin.close()
		if stdout:
			p.stdout.close()
		time.sleep(1)
		if p.wait() is None:
			print 'Could not terminate subprocess... forcing termination'
			p.kill()
			time.sleep(1)
			if p.wait() is None:
				print 'Could not kill subprocess... quitting'
				exit(os.EX_SOFTWARE)
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print 'The subprocess is still running. Killing it with os.kill'
			os.kill(pid, 15)
			try:
				os.kill(pid, 0)
			except OSError:
				pass
			else:
				print 'The process does not die... quitting'
				exit(os.EX_SOFTWARE)
		del pid
		p = None

def main():
	parser = OptionParser()
	parser.add_option('-r', '--runlist', dest='runlist', type='string', help='File containing a list of the RunSettings for each run')
	parser.add_option('-n', '--numcores', dest='numcores', type='int', default=2, help='number of runs to execute in parallel')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true', help='force to use the specified number of cores')

	(options, args) = parser.parse_args()
	runlist = options.runlist
	num = options.numcores
	force = options.force

	pp = RunNormalAnalysisParallel(runlist=runlist, num_cores=num, force=force)


if __name__ == '__main__':
	main()
