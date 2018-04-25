#!/usr/bin/env python
# import ROOT as ro
# import ipdb  # set_trace, launch_ipdb_on_exception
# from copy import deepcopy
import os, sys, shutil
sys.path.append('/afs/cern.ch/user/d/dsanzbec/rd42AnalysisTools')
from Utils import *
import subprocess as subp
import multiprocessing as mp
from ConfigParser import ConfigParser

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

# fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12}
ped_axis = {'min': 0, 'max': 2**12}
cm_axis = {'min': -100, 'max': 100}

class RD42Analysis:
	def __init__(self):
		print 'Creating RD42Analysis'
		self.run = 0
		self.total_events = 0
		self.dia_input = 0
		self.data_dir = ''
		self.out_dir = ''
		self.settings_dir = ''
		self.run_lists_dir = ''
		self.subdir = 'no_mask'
		self.do_even = False
		self.do_odd = False
		self.do_chs = False
		self.batch = True
		self.StripTelescopeAnalysis_path = '/afs/cern.ch/user/d/dsanzbec/StripTelescopeAnalysis'
		self.scratch_path = '/eos/user/d/dsanzbec/scratch/output'  # at lxplus
		# self.scratch_path = '/scratch/strip_telescope_tests/runDiego/output'  # at snickers

		self.first_event = 0
		self.num_events = 0
		self.do_pedestal = False
		self.do_cluster = False
		self.do_selection = False
		self.do_alignment = False
		self.do_transparent = False
		self.do_3d = False

		self.sub_pro = None

	def ReadInputFile(self, in_file=''):
		if in_file != '':
			if os.path.isfile(in_file):
				pars = ConfigParser()
				pars.read(in_file)
				print 'Loading job description from file:', in_file

				if pars.has_section('RUN'):
					if pars.has_option('RUN', 'StripTelescopeAnalysis_path'):
						self.StripTelescopeAnalysis_path = pars.get('RUN', 'StripTelescopeAnalysis_path')
					if pars.has_option('RUN', 'run'):
						self.run = pars.getint('RUN', 'run')
					else:
						ExitMessage('Must specify run under [RUN]. Exiting...')
					if pars.has_option('RUN', 'events'):
						self.total_events = pars.getint('RUN', 'events')
					else:
						ExitMessage('Must specify events under [RUN]. Exiting...')
					if pars.has_option('RUN', 'dia_input'):
						self.dia_input = pars.getint('RUN', 'dia_input')
					if pars.has_option('RUN', 'datadir'):
						self.data_dir = pars.get('RUN', 'datadir')
					else:
						ExitMessage('Must specify datadir under [RUN]. Exiting...')
					if pars.has_option('RUN', 'outputdir'):
						self.out_dir = pars.get('RUN', 'outputdir')
					else:
						ExitMessage('Must specify outputdir under [RUN]. Exiting...')
					if pars.has_option('RUN', 'settingsdir'):
						self.settings_dir = pars.get('RUN', 'settingsdir')
					else:
						ExitMessage('Must specify settingsdir under [RUN]. Exiting...')
					if pars.has_option('RUN', 'runlistsdir'):
						self.run_lists_dir = pars.get('RUN', 'runlistsdir')
					else:
						ExitMessage('Must specify runlistsdir under [RUN]. Exiting...')
					if pars.has_option('RUN', 'subdir'):
						self.subdir = pars.get('RUN', 'subdir')
					if pars.has_option('RUN', 'do_even'):
						self.do_even = pars.getboolean('RUN', 'do_even')
					if pars.has_option('RUN', 'do_odd'):
						self.do_odd = pars.getboolean('RUN', 'do_odd')
					if pars.has_option('RUN', 'do_chs'):
						self.do_chs = pars.getboolean('RUN', 'do_chs')
					if pars.has_option('RUN', 'batch'):
						self.batch = pars.getboolean('RUN', 'batch')

				if pars.has_section('ANALYSIS'):
					if pars.has_option('ANALYSIS', 'first_event'):
						self.first_event = pars.getint('ANALYSIS', 'first_event')
					if pars.has_option('ANALYSIS', 'num_events'):
						self.num_events = pars.getint('ANALYSIS', 'num_events')
					if pars.has_option('ANALYSIS', 'do_pedestal'):
						self.do_pedestal = pars.getboolean('ANALYSIS', 'do_pedestal')
					if pars.has_option('ANALYSIS', 'do_cluster'):
						self.do_cluster = pars.getboolean('ANALYSIS', 'do_cluster')
					if pars.has_option('ANALYSIS', 'do_selection'):
						self.do_selection = pars.getboolean('ANALYSIS', 'do_selection')
					if pars.has_option('ANALYSIS', 'do_alignment'):
						self.do_alignment = pars.getboolean('ANALYSIS', 'do_alignment')
					if pars.has_option('ANALYSIS', 'do_transparent'):
						self.do_transparent = pars.getboolean('ANALYSIS', 'do_transparent')
					if pars.has_option('ANALYSIS', 'do_3d'):
						self.do_3d = pars.getboolean('ANALYSIS', 'do_3d')

				self.num_events = self.total_events if self.num_events == 0 else self.num_events
				return
		ExitMessage('Input file "{i}" does not exist. Must input a valid file. Exiting'.format(i=in_file))

	def Create_Run_List(self):
		if not os.path.isdir(self.run_lists_dir):
			os.makedirs(self.run_lists_dir)
		ped = 1 if self.do_pedestal else 0
		clu = 1 if self.do_cluster else 0
		sele = 1 if self.do_selection else 0
		alig = 1 if self.do_alignment else 0
		tran = 1 if self.do_transparent else 0
		if not os.path.isdir(self.run_lists_dir+'/{f}'.format(f=self.subdir)):
			os.makedirs(self.run_lists_dir+'/{f}'.format(f=self.subdir))
		with open(self.run_lists_dir + '/{f}/RunList_'.format(f=self.subdir)+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t{al}\t0\t{t}\n#\n'.format(r=self.run, n=self.num_events, p=ped, c=clu, s=sele, al=alig, t=tran))
		if self.do_even or self.do_odd:
			if not os.path.isdir(self.run_lists_dir+'/'+self.subdir+'/odd'):
				os.makedirs(self.run_lists_dir+'/'+self.subdir+'/odd')
			with open(self.run_lists_dir+'/'+self.subdir+'/odd' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.num_events, p=ped, c=clu, s=sele))
			if not os.path.isdir(self.run_lists_dir+'/'+self.subdir+'/even'):
				os.makedirs(self.run_lists_dir+'/'+self.subdir+'/even')
			with open(self.run_lists_dir+'/'+self.subdir+'/even' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.num_events, p=ped, c=clu, s=sele))

	def Check_settings_file(self):
		if not os.path.isdir(self.settings_dir + '/' + self.subdir):
			os.makedirs(self.settings_dir + '/' + self.subdir)
		if not os.path.isfile(self.settings_dir + '/' + self.subdir + '/settings.{r}.ini'.format(r=self.run)):
			CreateDefaultSettingsFile(self.settings_dir + '/' + self.subdir, self.run, self.total_events, ev_ini=self.first_event, num_evs_ana=self.num_events, dia_input=self.dia_input)
		if self.do_even or self.do_odd:
			self.Copy_settings_to_even_odd()

	def Copy_settings_to_even_odd(self):
		if not os.path.isdir(self.settings_dir + '/'+self.subdir+'/even'):
			os.makedirs(self.settings_dir + '/'+self.subdir+'/even')
		if not os.path.isdir(self.settings_dir + '/'+self.subdir+'/odd'):
			os.makedirs(self.settings_dir + '/'+self.subdir+'/odd')
		shutil.copy(self.settings_dir + '/'+self.subdir+'/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/'+self.subdir+'/even/')
		shutil.copy(self.settings_dir + '/'+self.subdir+'/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/'+self.subdir+'/odd/')

	def CheckStripTelescopeAnalysis(self):
		if os.path.isdir(self.StripTelescopeAnalysis_path):
			if not os.path.isfile(self.StripTelescopeAnalysis_path + '/diamondAnalysis'):
				ExitMessage('{p}/diamondAnalysis does not exist. Exiting'.format(p=self.StripTelescopeAnalysis_path))
		else:
			ExitMessage('{d} does not exist. Exiting'.format(d=self.StripTelescopeAnalysis_path))

	def Print_subprocess_command(self, runlist, setting, outdir, inputdir):
		print 'Executing:\n{p}/diamondAnalysis -r {r} -s {s} -o {o} -i {i}\n'.format(p=self.StripTelescopeAnalysis_path, r=runlist, s=setting, o=outdir, i=inputdir)

	def RunNormalAnalysis(self):
		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.run))
		RecreateSoftLink(self.out_dir + '/' + self.subdir + '/' + str(self.run), self.scratch_path, str(self.run) + '_' + self.subdir)
		self.Print_subprocess_command('{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), self.settings_dir + '/' + self.subdir, self.out_dir + '/' + self.subdir, self.data_dir + '/' + str(self.run))
		if self.batch:
			self.sub_pro = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), close_fds=True)
		else:
			self.sub_pro = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		while self.sub_pro.poll() is None:
			pass
		if self.sub_pro.poll() == 0:
			print 'Run finished successfully'
		else:
			print 'Run could have failed. Obtained return code:', self.sub_pro.poll()
		CloseSubprocess(self.sub_pro, True, False)
		if self.do_odd:
			CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/odd/' + str(self.run))
			self.Print_subprocess_command('{d}/{sd}/odd/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), self.settings_dir + '/' + self.subdir + '/odd', self.out_dir + '/' + self.subdir + '/odd', self.data_dir + '/' + str(self.run))
			self.sub_pro = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/odd/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir + '/odd', '-o', self.out_dir + '/' + self.subdir + '/odd', '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), close_fds=True)
			while self.sub_pro.poll() is None:
				pass
			if self.sub_pro.poll() == 0:
				print 'Run odd finished'
			else:
				print 'Run odd could have failed. Obtained return code:', self.sub_pro.poll()
			CloseSubprocess(self.sub_pro, True, False)
		if self.do_even:
			CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/even/' + str(self.run))
			self.Print_subprocess_command('{d}/{sd}/even/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), self.settings_dir + '/' + self.subdir + '/even', self.out_dir + '/' + self.subdir + '/even', self.data_dir + '/' + str(self.run))
			self.sub_pro = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/even/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir + '/even', '-o', self.out_dir + '/' + self.subdir + '/even', '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), close_fds=True)
			while self.sub_pro.poll() is None:
				pass
			if self.sub_pro.poll() == 0:
				print 'Run even finished'
			else:
				print 'Run even could have failed. Obtained return code:', self.sub_pro.poll()
			CloseSubprocess(self.sub_pro, True, False)


def main(argv):
	print 'arguments are:', argv
	infile = argv[1]
	rd42 = RD42Analysis()
	rd42.ReadInputFile(infile)
	rd42.Create_Run_List()
	rd42.Check_settings_file()
	rd42.CheckStripTelescopeAnalysis()
	rd42.RunNormalAnalysis()
	print '\nFinished :)\n'
	sys.exit(os.EX_OK)

if __name__ == '__main__':
	main(sys.argv)
