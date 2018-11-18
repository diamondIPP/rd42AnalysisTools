#!/usr/bin/env python
from optparse import OptionParser
import os, sys, shutil
from Utils import *
import glob
import ipdb

__author__ = 'DA'


class CreateRunSettings:
	def __init__(self, dir='RunSettings', run=0, softwaredir='/sdvlp', inputdir='~/data'):
		print 'Creating run settings file for run', run
		self.run = run
		self.dir = os.path.abspath(os.path.expanduser(os.path.expandvars(dir)))
		self.softwaredir = os.path.abspath(os.path.expanduser(os.path.expandvars(softwaredir)))
		self.STAdir = self.softwaredir
		self.inputdir = os.path.abspath(os.path.expanduser(os.path.expandvars(inputdir)))
		self.ctcpath = self.softwaredir
		self.events = 0
		self.diainput = 1
		self.dutname = 'diamond'
		self.voltage = 0
		self.saturation = 4095
		self.eta_corr_limit = 0.0005
		self.max_transp_clust_size = 10
		self.num_strips_transp = 5
		self.chi2 = 5
		self.transparentChi2 = 5
		self.transparentAlignment = 0
		self.ph_dia_max = 4000
		self.ph_num_bins = 200

		CreateDirectoryIfNecessary(self.dir)

		self.STAdir = glob.glob('{d}/StripTelescopeAnalysis'.format(d=self.softwaredir))
		self.STAdir = self.ValidatePath(self.STAdir, 'Enter the full path to the StripTelescopeAnalysis directory: ', isdir=True)

		self.ctcpath = glob.glob('{d}/rd42AnalysisTools/CrossTalkCorrections/FeedThroughCorrection'.format(d=self.softwaredir))
		self.ctcpath = self.ValidatePath(self.ctcpath, 'Enter the full path to the FeedThroughCorrection executable: ', isdir=False)

		self.run = self.ValidateInteger(self.run, 'run = {r}.\tPress ENTER if correct, otherwise enter a new run number: '.format(r=self.run))
		self.events = self.ValidateInteger(self.events, 'Total events = {ev}.\tPress ENTER if correct, otherwise enter a new number for total events in the run: '.format(ev=self.events))
		self.diainput = self.ValidateInteger(self.diainput, 'diamond input = {d}.\tPress ENTER if correct, otherwise enter 0 if it was in sirrocco 4 or 1 if it was in sirrocco 5: '.format(d=self.diainput), [0, 1])
		self.dutname = self.ValidateString(self.dutname, 'dut name = {n}.\tPress ENTER if correct, otherwise enter a new dut name: '.format(n=self.dutname))
		self.voltage = self.ValidateFloat(self.voltage, 'Bias voltage = {v}.\tPress ENTER if correct, otherwise enter a new bias voltage: '.format(v=self.voltage))
		self.saturation = 4095 if self.diainput == 1 else 3367
		self.eta_corr_limit = self.ValidateFloat(self.eta_corr_limit, 'eta correction convergence limit = {v}.\tPress ENTER if correct, otherwise enter a eta correction convergence limit: '.format(v=self.eta_corr_limit))
		self.max_transp_clust_size = self.ValidateInteger(self.max_transp_clust_size, 'Transparent cluster size = {v}.\tPress ENTER if correct, otherwise enter a new transparent cluster size: '.format(v=self.max_transp_clust_size))
		self.num_strips_transp = self.ValidateInteger(self.num_strips_transp, 'Number of strips within the transparent cluster for charge = {v}.\tPress ENTER if correct, otherwise enter the number of strips to calculate the charge in the transparent cluster: '.format(v=self.num_strips_transp))
		self.chi2 = self.ValidateFloat(self.chi2, 'chi2 = {v}.\tPress ENTER if correct, otherwise enter the maximum chi2 value to consider: '.format(v=self.chi2))
		self.transparentChi2 = self.ValidateFloat(self.transparentChi2, 'transparentChi2 = {v}.\tPress ENTER if correct, otherwise enter the maximum transparent chi2 value to consider: '.format(v=self.transparentChi2))
		self.transparentAlignment = self.ValidateInteger(self.transparentAlignment, 'Do transparent alignment = {v}.\tPress ENTER if correct, otherwise enter 0 or 1 depending if you don\'t want or do want to do transparent alignment: '.format(v=bool(self.transparentAlignment)))
		self.ph_dia_max = self.ValidateFloat(self.ph_dia_max, 'Max ph value = {v}.\tPress ENTER if correct, otherwise enter a new max for ph plots: '.format(v=self.ph_dia_max))
		self.ph_num_bins = self.ValidateInteger(self.ph_num_bins, 'Num bins for ph = {v}.\tPress ENTER if correct, otherwise enter a new number of bins for ph plots: '.format(v=self.ph_num_bins))

		self.inputrundir = glob.glob('{d}/*/cern_RD42_*/{r}'.format(d=self.inputdir, r=self.run))
		self.inputrundir = self.ValidatePath(self.inputrundir, 'Enter the full path to the run directory with the raw files are located: ', True)
		self.inputrundir = self.inputrundir.split('/{r}'.format(r=self.run))[0]

		self.outputdir = self.inputrundir.split('/cern_RD42')[0] + '/output'
		self.outputdir = self.ValidatePath([self.outputdir], 'output dir = {d}.\tPress ENTER if correct, otherwise enter a new full path for the output directory: ', couldbenew=True)

		self.settingsdir = os.path.abspath(os.path.expanduser('~/settings'))
		self.settingsdir = self.ValidatePath([self.settingsdir], 'settings path = {v}.\tPress ENTER if correct, otherwise enter a new path for the analysis settings files: ', couldbenew=True)

		self.runlistdir = os.path.abspath(os.path.expanduser('~/RunLists'))
		self.runlistdir = self.ValidatePath([self.runlistdir], 'RunLists path = {v}.\tPress ENTER if correct, otherwise enter a new path for the analysis RunLists files: ', couldbenew=True)

		self.scratch = os.path.abspath(os.path.expanduser('~/scratch'))

		self.subdir = 'no_mask'
		self.subdir = self.ValidateString(self.subdir, 'subdir = {v}.\tPress ENTER if correct, otherwise enter a new subdir string: '.format(v=self.subdir))

		self.filename = 'settings{r}_{sd}'.format(r=self.run, sd=self.subdir)
		self.filename = self.ValidateString(self.filename, 'run settings file name = {f}.ini.\tPress ENTER if correct, otherwise enter a new filenme for the run settings file: '.format(f=self.filename))

	def ValidatePath(self, path, message, isdir=True, couldbenew=False):
		while len(path) != 1:
			path = [raw_input(message)]
			if isdir:
				if not os.path.isdir(path[0]) and not couldbenew:
					path = []
			else:
				if not os.path.isfile(path[0]) and not couldbenew:
					path = []
		return path[0]

	def ValidateInteger(self, var, message, options=[]):
		cont = False
		while not cont:
			temp = raw_input(message)
			if temp == '':
				if len(options) != 0:
					if var in options:
						val = var
						cont = True
				else:
					val = var
					cont = True
			else:
				if IsInt(temp):
					val = int(temp)
					if len(options) != 0:
						if val in options:
							cont = True
					else:
						cont = True
		return val

	def ValidateFloat(self, var, message, options=[]):
		cont = False
		while not cont:
			temp = raw_input(message)
			if temp == '':
				if len(options) != 0:
					if var in options:
						val = var
						cont = True
				else:
					val = var
					cont = True
			else:
				if IsFloat(temp):
					val = int(temp)
					if len(options) != 0:
						if val in options:
							cont = True
					else:
						cont = True
		return val

	def ValidateString(self, var, message):
		temp = raw_input(message)
		val = var
		if temp != '':
			val = temp
		return val

	def CreateRunSettingsFile(self):
		with open('{d}/{f}.ini'.format(d=self.dir, f=self.filename), 'w') as f0:
			f0.write('[RUN]\n')
			f0.write('StripTelescopeAnalysis_path = {v}\n'.format(v=self.STAdir))
			f0.write('cross_talk_correction_path = {v}\n'.format(v=self.ctcpath))
			f0.write('run = {v}\n'.format(v=self.run))
			f0.write('events = {v}\n'.format(v=self.events))
			f0.write('dia_input = {v}\n'.format(v=self.diainput))
			f0.write('dut_name = {v}\n'.format(v=self.dutname))
			f0.write('dut_volt = {v}\n'.format(v=self.voltage))
			f0.write('dia_saturation = {v}\n'.format(v=self.saturation))
			f0.write('eta_corr_limit = {v}\n'.format(v=self.eta_corr_limit))
			f0.write('max_transparent_cluster_size = {v}\n'.format(v=self.max_transp_clust_size))
			f0.write('num_highest_transparent_cluster = {v}\n'.format(v=self.num_strips_transp))
			f0.write('chi2 = {v}\n'.format(v=self.chi2))
			f0.write('transparentChi2 = {v}\n'.format(v=self.transparentChi2))
			f0.write('transparentAlignment = {v}\n'.format(v=int(self.transparentAlignment)))
			f0.write('ph_dia_max = {v}\n'.format(v=self.ph_dia_max))
			f0.write('ph_num_bins = {v}\n'.format(v=self.ph_num_bins))
			f0.write('datadir = {v}\n'.format(v=self.inputrundir))
			f0.write('outputdir = {v}\n'.format(v=self.outputdir))
			f0.write('settingsdir = {v}\n'.format(v=self.settingsdir))
			f0.write('runlistsdir = {v}\n'.format(v=self.runlistdir))
			f0.write('scratch_path = {v}\n'.format(v=self.scratch))
			f0.write('subdir = {v}\n'.format(v=self.subdir))
			f0.write('do_even = {v}\n'.format(v='False'))
			f0.write('do_odd = {v}\n'.format(v='False'))
			f0.write('do_chs = {v}\n'.format(v='False'))
			f0.write('batch = {v}\n'.format(v='False'))
			f0.write('symlinks = {v}\n'.format(v='True'))
			f0.write('delete_old = {v}\n\n'.format(v='False'))
			f0.write('[ANALYSIS]\n')
			f0.write('first_event = {v}\n'.format(v=0))
			f0.write('num_events = {v}\n'.format(v=0))
			f0.write('do_pedestal = {v}\n'.format(v=0))
			f0.write('do_cluster = {v}\n'.format(v=0))
			f0.write('do_selection = {v}\n'.format(v=0))
			f0.write('do_alignment = {v}\n'.format(v=0))
			f0.write('do_transparent = {v}\n'.format(v=0))
			f0.write('do_3d = {v}\n'.format(v=0))
			f0.write('do_cross_talk_calc = {v}\n'.format(v=0))
		print 'Finished creating the run settings file {f}.ini under {d} path :)'.format(f=self.filename, d=self.dir)
		exit(os.EX_OK)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', default=0, type='int', help='Run to be analysed (e.g. 22022)')
	parser.add_option('-d', '--dir', dest='dir', default='RunSettings', type='string', help='Directory where the run settings file will be saved (e.g. RunSettings)')
	parser.add_option('-s', '--softwaredir', dest='softwaredir', default='/sdvlp', type='string', help='Path to directory where the different analysis softwares are located (e.g. /sdvlp)')
	parser.add_option('-i', '--inputdir', dest='inputdir', default='~/data', type='string', help='Path to directory where the years of the testbeam campaings containing the raw files are located (e.g. ~/data)')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	softwaredir = str(options.softwaredir)
	inputdir = str(options.inputdir)

	crs = CreateRunSettings(dir=dir, run=run, softwaredir=softwaredir, inputdir=inputdir)
	crs.CreateRunSettingsFile()
