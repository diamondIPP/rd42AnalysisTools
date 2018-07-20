#!/usr/bin/env python
import time
import subprocess as subp
import multiprocessing as mp
import os, shutil
import ipdb
from optparse import OptionParser

class RunNormalAnalysisParallel:
	def __init__(self, runlist, num_cores=2):
		self.runlist = runlist
		self.num_cores = num_cores if num_cores <= int(mp.cpu_count()/2.0) else int(mp.cpu_count()/2.0)
		self.num_runs = 0
		self.settings_list = []
		self.job_chunks = []
		self.analysis_processes = {}
		self.workind_dir = os.getcwd()
		self.queue = {}
		self.queue_running = {}
		if not os.path.isfile(self.runlist):
			print 'File', self.runlist, 'does not exist. Exiting'
			exit(os.EX_CONFIG)
		print 'Starting parallel analysis using runlist', self.runlist, 'and using', self.num_cores, 'cores simultaneously'
		self.ReadRunList()
		self.RunParallelAnalysis()

	def ReadRunList(self):
		with open(self.runlist, 'r') as rl:
			lines = rl.readlines()
			self.settings_list = [line.split('\n')[0] for line in lines if os.path.isfile(line.split('\n')[0])]
			self.num_runs = len(self.settings_list)
			self.job_chunks = [self.settings_list[i:i + self.num_cores] for i in xrange(0, self.num_runs, self.num_cores)]
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
		self.queue = {i: None for i in xrange(self.num_cores)}
		self.queue_running = {i: False for i in xrange(self.num_cores)}
		with open(os.devnull, 'w') as FNULL:
			do_add_queue = True
			for jobi in self.settings_list:
				if do_add_queue:
					pos_q, nfree = self.GetAvailablePos()
					if nfree > 0:
						if nfree == 1:
							print 'Showing output for run', jobi
							self.queue[pos_q] = subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(jobi), '--normal'], bufsize=-1, stdin=subp.PIPE, close_fds=True)
							do_add_queue = False
						else:
							self.queue[pos_q] = subp.Popen(['{wd}/rd42AnalysisBatch.py'.format(wd=self.workind_dir), '-w', os.getcwd(), '-s', os.path.abspath(jobi), '--normal'], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True)
						self.queue_running[pos_q] = True
					else:
						do_add_queue = False
				# if not do_add_queue:
				# 	if np.all


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

	(options, args) = parser.parse_args()
	runlist = options.runlist
	num = options.numcores

	pp = RunNormalAnalysisParallel(runlist=runlist, num_cores=num)


if __name__ == '__main__':
	main()
