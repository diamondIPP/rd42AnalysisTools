#!/usr/bin/env python
from optparse import OptionParser
import progressbar
import os, sys, shutil
from Utils import *

__author__ = 'DA'

option_list = ['raw', 'pedestal', 'cluster', 'selection', 'alignment', 'transparent']

link_list = {'raw': ['rawData.{r}.root'], 'pedestal': ['pedestalData.{r}.root', 'pedestalAnalysis'], 'cluster': ['clusterData.{r}.root', 'clustering', 'etaCorrection.{r}.root'], 'selection': ['selectionData.{r}.root', 'selections', 'selectionAnalysis'], 'alignment': ['alignment.{r}.root', 'alignment'], 'transparent': ['transparentAnalysis']}
copy_list = {'cluster': ['crossTalkCorrectionFactors.{r}.txt']}

class LinkFiles:
	def __init__(self, source_subdir='', dest_subdir='', run=0, upto='alignment', force=False):
		print 'Linking files for run:', run, 'upto', upto
		self.run = run
		self.upto = upto
		self.s_subdir = source_subdir
		self.d_subdir = dest_subdir
		self.force = force
		self.bar = None
		self.cumulative = [option_list[:i+1] for i in xrange(len(option_list)) if option_list[i] == self.upto][0]

		if not os.path.isdir(self.d_subdir + '/' + str(self.run)):
			os.makedirs(self.d_subdir + '/' + str(self.run))

		os.chdir(self.d_subdir + '/' + str(self.run))

		self.relative_path = os.path.relpath(self.s_subdir + '/' + str(self.run))

		self.CreateProgressBar(len(self.cumulative))
		if not self.bar:
			self.bar.start()
		count = 0
		for keyc in self.cumulative:
			for elem in link_list[keyc]:
				do_link = False
				if os.path.isdir(self.relative_path + '/' + elem.format(r=self.run)) or os.path.isfile(self.relative_path + '/' + elem.format(r=self.run)) or os.path.islink(self.relative_path + '/' + elem.format(r=self.run)):
					if os.path.isdir(elem.format(r=self.run)) or os.path.isfile(elem.format(r=self.run)) or os.path.islink(elem.format(r=self.run)):
						do_link = False
						if self.force:
							if os.path.isdir(elem.format(r=self.run)):
								shutil.rmtree(elem.format(r=self.run), True)
							else:
								os.unlink(elem.format(r=self.run))
							do_link = True
					else:
						do_link = True
					if do_link:
						os.symlink(self.relative_path + '/' + elem.format(r=self.run), elem.format(r=self.run))
			if keyc in copy_list.keys():
				for elem in copy_list[keyc]:
					if os.path.isdir(self.relative_path + '/' + elem.format(r=self.run)) or os.path.isfile(self.relative_path + '/' + elem.format(r=self.run)) or os.path.islink(self.relative_path + '/' + elem.format(r=self.run)):
						if os.path.isdir(elem.format(r=self.run)) or os.path.isfile(elem.format(r=self.run)) or os.path.islink(elem.format(r=self.run)):
							do_copy = False
							if self.force:
								if os.path.isdir(elem.format(r=self.run)):
									shutil.rmtree(elem.format(r=self.run), True)
								else:
									os.unlink(elem.format(r=self.run))
								do_copy = True
						else:
							do_copy = True
						if do_copy:
							if os.path.isdir(self.relative_path + '/' + elem.format(r=self.run)):
								shutil.copytree(self.relative_path + '/' + elem.format(r=self.run), elem.format(r=self.run), True)
							else:
								shutil.copy2(self.relative_path + '/' + elem.format(r=self.run), elem.format(r=self.run))
			count += 1
			if not self.bar:
				self.bar.update(count)
		if not self.bar:
			self.bar.finish()

		print 'Finished linking files upto', self.upto, 'for run', self.run, ':)'
		exit(os.EX_OK)

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', default=22022, type='int', help='Run to be analysed (e.g. 22022)')
	parser.add_option('-s', '--sourcesubdir', dest='sdir', default='.', type='string', help='Source subdir containing the run folder')
	parser.add_option('-d', '--destsubdir', dest='ddir', default='.', type='string', help='Destination subdir where the run folder will have links')
	parser.add_option('-u', '--upto', dest='upto', default='alignment', type='string', help='The input can be "raw", "pedestal", "cluster", "selection", "alignment", "transparent". Example: If "cluster" is entered, then it will only link "raw", "pedestal" and "cluster"')

	(options, args) = parser.parse_args()
	run = int(options.run)
	sdir = str(options.sdir)
	ddir = str(options.ddir)
	upto = str(options.upto)
	if upto not in option_list:
		print '-u (--upto) option must be one of these: "raw", "pedestal", "cluster", "selection", "alignment", "transparent". Exiting'
		exit(os.EX_CONFIG)

	if not os.path.isdir(sdir):
		print 'Source subdirectory does not exist. Exiting'
		exit(os.EX_CONFIG)

	if not os.path.isdir(sdir + '/' + str(run)):
		print 'Source subdirectory does not have run folder', run, 'exiting.'
		exit(os.EX_CONFIG)

	sdir = os.path.abspath(sdir)
	ddir = os.path.abspath(ddir)

	lf = LinkFiles(source_subdir=sdir, dest_subdir=ddir, run=run, upto=upto)
