#!/usr/bin/env python
from optparse import OptionParser
import progressbar
import os, sys, shutil
from Utils import *

__author__ = 'DA'

option_list = ['raw', 'pedestal', 'cluster', 'selection', 'alignment', 'transparent']

list1 = {'raw': ['rawData.{r}.root'], 'pedestal': ['pedestalData.{r}.root', 'pedestalAnalysis'], 'cluster': ['clusterData.{r}.root', 'clustering', 'etaCorrection.{r}.root'], 'selection': ['selectionData.{r}.root', 'selections', 'selectionAnalysis'], 'alignment': ['alignment.{r}.root', 'alignment'], 'transparent': ['transparentAnalysis', 'transparent.{r}.root','transparent2.{r}.root']}
list2 = {'cluster': ['crossTalkCorrectionFactors.{r}.txt']}

class UnlinkFiles:
	def __init__(self, source_subdir='', run=0, element_list=[], force=False):
		print 'Uninking/deleting files and folders of the categories', element_list, 'for run:', run
		self.run = run
		self.element_list = element_list
		self.s_subdir = source_subdir
		self.force = force
		self.bar = None

		self.CreateProgressBar(len(self.element_list))
		if not self.bar:
			self.bar.start()
		count = 0
		for keyc in self.element_list:
			for elem in list1[keyc]:
				path_elem = self.s_subdir + '/' + str(self.run) + '/' + elem.format(r=self.run)
				if os.path.islink(path_elem):
					os.unlink(path_elem)
				elif os.path.isdir(path_elem):
					shutil.rmtree(path_elem)
				elif os.path.isfile(path_elem):
					os.remove(path_elem)
			if keyc in list2.keys():
				for elem in list2[keyc]:
					path_elem = self.s_subdir + '/' + str(self.run) + '/' + elem.format(r=self.run)
					if os.path.islink(path_elem):
						os.unlink(path_elem)
					elif os.path.isdir(path_elem):
						shutil.rmtree(path_elem)
					elif os.path.isfile(path_elem):
						os.remove(path_elem)
			count += 1
			if not self.bar:
				self.bar.update(count)
		if not self.bar:
			self.bar.finish()

		print 'Finished unlinking/deleting files of the caterogires', self.element_list, 'for run', self.run, ':)'
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
	parser.add_option('-e', '--elements', dest='elem', default='[]', type='string', help='List of elements contained in square brackets of the files and directories that will be unlinked/deleted. The options are: "raw", "pedestal", "cluster", "selection", "alignment", "transparent". e.g. [cluster,transparent,raw]')

	(options, args) = parser.parse_args()
	run = int(options.run)
	sdir = str(options.sdir)
	elemtemp = str(options.elem)
	elem = elemtemp.split('[')[1].split(']')[0].split(',')
	if not set(elem).issubset(option_list):
		print '-u (--upto) option must be one of these: "raw", "pedestal", "cluster", "selection", "alignment", "transparent". Exiting'
		exit(os.EX_CONFIG)

	if not os.path.isdir(sdir):
		print 'Source subdirectory does not exist. Exiting'
		exit(os.EX_CONFIG)

	if not os.path.isdir(sdir + '/' + str(run)):
		print 'Source subdirectory does not have run folder', run, 'exiting.'
		exit(os.EX_CONFIG)

	sdir = os.path.abspath(sdir)

	lf = UnlinkFiles(source_subdir=sdir, run=run, element_list=elem)
