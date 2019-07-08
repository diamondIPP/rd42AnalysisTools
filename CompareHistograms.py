#!/usr/bin/env python
from optparse import OptionParser
import ROOT as ro
import progressbar
import os, sys, shutil
from Utils import *

__author__ = 'DA'


class CompareHistograms:
	def __init__(self, source_subdir='', dest_subdir='', pth='', name=''):
		print 'Comparing at file:', name, 'inside', source_subdir, 'and', dest_subdir
		self.name = name
		self.pth = pth
		self.s_subdir = source_subdir
		self.d_subdir = dest_subdir
		self.bar = None

		if not os.path.isdir(self.d_subdir + '/' + self.pth):
			ExitMessage(self.d_subdir + '/' + self.pth + ' does not exist', os.EX_DATAERR)

		if not os.path.isdir(self.s_subdir + '/' + self.pth):
			ExitMessage(self.s_subdir + '/' + self.pth + ' does not exist', os.EX_DATAERR)

		if not os.path.isfile(self.d_subdir + '/' + self.pth + '/' + self.name + '.root'):
			ExitMessage(self.d_subdir + '/' + self.pth + '/' + self.name + '.root does not exist', os.EX_DATAERR)

		if not os.path.isfile(self.s_subdir + '/' + self.pth + '/' + self.name + '.root'):
			ExitMessage(self.s_subdir + '/' + self.pth + '/' + self.name + '.root does not exist', os.EX_DATAERR)

		self.fs = ro.TFile(self.s_subdir + '/' + self.pth + '/' + self.name + '.root', 'read')
		self.fd = ro.TFile(self.d_subdir + '/' + self.pth + '/' + self.name + '.root', 'read')



if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-s', '--sourcedir', dest='sdir', default='.', type='string', help='Source subdir containing the run folder')
	parser.add_option('-d', '--destdir', dest='ddir', default='.', type='string', help='Destination subdir where the run folder will have links')
	parser.add_option('-p', '--pathinternal', dest='path', default='', type='string', help='Internal path to the root file for both runs')
	parser.add_option('-n', '--name', dest='name', default='', type='string', help='Name of the root file to compare (without .root extension)')

	(options, args) = parser.parse_args()
	sdir = str(options.sdir)
	ddir = str(options.ddir)
	pth = str(options.path)
	name = str(options.name)

	sdir = os.path.abspath(sdir)
	ddir = os.path.abspath(ddir)

	ch = CompareHistograms(source_subdir=sdir, dest_subdir=ddir, pth=pth, name=name)
