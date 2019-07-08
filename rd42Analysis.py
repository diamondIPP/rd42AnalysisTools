#!/usr/bin/env python
import os, sys, shutil
sys.path.append('/home/sandiego/rd42AnalysisTools')  # TODO: HARDCODED!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from ConfigParser import ConfigParser
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
# from NoiseExtraction import NoiseExtraction
from Utils import *
import subprocess as subp
import multiprocessing as mp
import time

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12}
ped_axis = {'min': 0, 'max': 2**12}
cm_axis = {'min': -100, 'max': 100}


class RD42Analysis:
	def __init__(self, working_dir='.', verb=False):
		print 'Creating RD42Analysis'
		self.run = 0
		self.total_events = 0
		self.run1, self.run2, self.total_events1, self.total_events2, self.merged_run_number, self.num_events_merged, self.do_merged = 0, 0, 0, 0, 0, 0, False
		self.dia_input = 0
		self.ped_buffer = 500
		self.dut_name = 'default'
		self.dut_volt = 0
		self.dia_saturation = 4095
		self.eta_corr_limit = 0.0005
		self.max_transparent_cluster_size = 10
		self.num_highest_transparent_cluster = 5
		self.chi2 = None
		self.trans_chi2 = None
		self.trans_align = False
		self.data_dir = ''
		self.out_dir = ''
		self.settings_dir = ''
		self.run_lists_dir = ''
		self.subdir = 'no_mask'
		self.do_even = False
		self.do_odd = False
		self.do_chs = False
		self.batch = not verb
		self.StripTelescopeAnalysis_path = '/afs/cern.ch/user/d/dsanzbec/StripTelescopeAnalysis'
		self.scratch_path = '/eos/user/d/dsanzbec/scratch/output'  # at lxplus
		self.cross_talk_correction_path = None
		self.symlinks = True
		self.ph_dia_max = 4096
		self.ph_dia_bins = 512
		self.yOffset = 0
		# self.scratch_path = '/scratch/strip_telescope_tests/runDiego/output'  # at snickers

		self.delete_old = False
		self.first_event = 0
		self.num_events = 0
		self.alignment_events = None
		self.do_pedestal = False
		self.do_cluster = False
		self.do_cross_talk = False
		self.do_selection = False
		self.do_alignment = False
		self.do_transparent = False
		self.do_3d = False
		self.current_dir = os.getcwd()
		self.working_dir = os.getcwd()

		self.sub_pro, self.sub_pro1, self.sub_pro2, self.sub_pro_crosstalk = None, None, None, None
		# self.process_f = None
		# self.process_e = None
		# self.process_o = None
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(999)

	def ReadInputFile(self, in_file=''):
		if in_file != '':
			if os.path.isfile(in_file):
				pars = ConfigParser()
				pars.read(in_file)
				print 'Loading job description from file:', in_file

				if pars.has_section('RUN'):
					if pars.has_option('RUN', 'StripTelescopeAnalysis_path'):
						self.StripTelescopeAnalysis_path = pars.get('RUN', 'StripTelescopeAnalysis_path')
					if pars.has_option('RUN', 'cross_talk_correction_path'):
						self.cross_talk_correction_path = pars.get('RUN', 'cross_talk_correction_path')
					if pars.has_option('RUN', 'run'):
						self.run = pars.getint('RUN', 'run')
					else:
						ExitMessage('Must specify run under [RUN]. Exiting...')
					if pars.has_option('RUN', 'run2') and pars.has_option('RUN', 'events2') and pars.has_option('RUN', 'merged_run_number'):
						self.run1 = pars.getint('RUN', 'run')
						self.run2 = pars.getint('RUN', 'run2')
						self.total_events1 = pars.getint('RUN', 'events')
						self.total_events2 = pars.getint('RUN', 'events2')
						self.merged_run_number = pars.getint('RUN', 'merged_run_number')
						self.do_merged = True
						if pars.has_section('ANALYSIS'):
							if pars.has_option('ANALYSIS', 'num_events'):
								self.num_events_merged = pars.getint('ANALYSIS', 'num_events')
					if pars.has_option('RUN', 'events'):
						self.total_events = pars.getint('RUN', 'events')
					else:
						ExitMessage('Must specify events under [RUN]. Exiting...')
					if pars.has_option('RUN', 'alignment_events'):
						self.alignment_events = pars.getint('RUN', 'alignment_events')
					if pars.has_option('RUN', 'dut_name'):
						self.dut_name = pars.get('RUN', 'dut_name')
					if pars.has_option('RUN', 'dut_volt'):
						self.dut_volt = pars.get('RUN', 'dut_volt')
					if pars.has_option('RUN', 'dia_input'):
						self.dia_input = pars.getint('RUN', 'dia_input')
					if pars.has_option('RUN', 'dia_saturation'):
						self.dia_saturation = pars.getint('RUN', 'dia_saturation')
					if pars.has_option('RUN', 'eta_corr_limit'):
						self.eta_corr_limit = pars.getfloat('RUN', 'eta_corr_limit')
					if pars.has_option('RUN', 'max_transparent_cluster_size'):
						self.max_transparent_cluster_size = pars.getint('RUN', 'max_transparent_cluster_size')
					if pars.has_option('RUN', 'num_highest_transparent_cluster'):
						self.num_highest_transparent_cluster = pars.getint('RUN', 'num_highest_transparent_cluster')
					if pars.has_option('RUN', 'chi2'):
						self.chi2 = pars.getfloat('RUN', 'chi2')
					if pars.has_option('RUN', 'transparentChi2'):
						self.trans_chi2 = pars.getfloat('RUN', 'transparentChi2')
					if pars.has_option('RUN', 'transparentAlignment'):
						self.trans_align = pars.getboolean('RUN', 'transparentAlignment')
					if pars.has_option('RUN', 'ph_dia_max'):
						self.ph_dia_max = pars.getint('RUN', 'ph_dia_max')
					if pars.has_option('RUN', 'ph_num_bins'):
						self.ph_dia_bins = pars.getint('RUN', 'ph_num_bins')
					if pars.has_option('RUN', 'yOffset'):
						self.yOffset = pars.getfloat('RUN', 'yOffset')
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
					if pars.has_option('RUN', 'symlinks'):
						self.symlinks = pars.getboolean('RUN', 'symlinks')
					if pars.has_option('RUN', 'delete_old'):
						self.delete_old = pars.getboolean('RUN', 'delete_old')
					if pars.has_option('RUN', 'ped_buffer'):
						self.ped_buffer = pars.getint('RUN', 'ped_buffer')

				if pars.has_section('ANALYSIS'):
					if pars.has_option('ANALYSIS', 'first_event'):
						self.first_event = pars.getint('ANALYSIS', 'first_event')
					if pars.has_option('ANALYSIS', 'num_events'):
						self.num_events = pars.getint('ANALYSIS', 'num_events')
					if pars.has_option('ANALYSIS', 'do_pedestal'):
						self.do_pedestal = pars.getboolean('ANALYSIS', 'do_pedestal')
					if pars.has_option('ANALYSIS', 'do_cluster'):
						self.do_cluster = pars.getboolean('ANALYSIS', 'do_cluster')
					if pars.has_option('ANALYSIS', 'do_cross_talk_calc'):
						if self.cross_talk_correction_path is None:
							print 'Will not do feed through (cross_talk) correction because the path to the executable has not been specified in cross_talk_correction_path.'
						else:
							self.do_cross_talk = pars.getboolean('ANALYSIS', 'do_cross_talk_calc')
					if pars.has_option('ANALYSIS', 'do_selection'):
						self.do_selection = pars.getboolean('ANALYSIS', 'do_selection')
					if pars.has_option('ANALYSIS', 'do_alignment'):
						self.do_alignment = pars.getboolean('ANALYSIS', 'do_alignment')
					if pars.has_option('ANALYSIS', 'do_transparent'):
						self.do_transparent = pars.getboolean('ANALYSIS', 'do_transparent')
					if pars.has_option('ANALYSIS', 'do_3d'):
						self.do_3d = pars.getboolean('ANALYSIS', 'do_3d')

				self.num_events = self.total_events if self.num_events == 0 else self.num_events
				self.alignment_events = RoundInt(self.num_events / 10.0) if not self.alignment_events else self.alignment_events
				return
		ExitMessage('Input file "{i}" does not exist. Must input a valid file. Exiting'.format(i=in_file))

	def Convert_Files(self):
		CreateDirectoryIfNecessary(self.run_lists_dir)
		ped1 = clu1 = sele1 = alig1 = tran1 = ped2 = clu2 = sele2 = alig2 = tran2 = 0
		CreateDirectoryIfNecessary(self.run_lists_dir+'/{f}'.format(f=self.subdir))
		with open(self.run_lists_dir + '/{f}/RunList_'.format(f=self.subdir)+str(self.run1)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t{al}\t0\t{t}\n#\n'.format(r=self.run1, n=self.total_events1, p=ped1, c=clu1, s=sele1, al=alig1, t=tran1))
		with open(self.run_lists_dir + '/{f}/RunList_'.format(f=self.subdir)+str(self.run2)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t{al}\t0\t{t}\n#\n'.format(r=self.run2, n=self.total_events2, p=ped2, c=clu2, s=sele2, al=alig2, t=tran2))
		CreateDirectoryIfNecessary(self.settings_dir + '/' + self.subdir)
		CreateDefaultSettingsFile(self.settings_dir + '/' + self.subdir, self.run1, self.total_events1, dut_name=self.dut_name, dut_volt=self.dut_volt, ev_ini=0, num_evs_ana=self.total_events1, dia_input=self.dia_input, dia_sat=self.dia_saturation, max_trans_clust=self.max_transparent_cluster_size, num_highest_trans=self.num_highest_transparent_cluster, chi2=self.chi2)
		CreateDefaultSettingsFile(self.settings_dir + '/' + self.subdir, self.run2, self.total_events2, dut_name=self.dut_name, dut_volt=self.dut_volt, ev_ini=0, num_evs_ana=self.total_events2, dia_input=self.dia_input, dia_sat=self.dia_saturation, max_trans_clust=self.max_transparent_cluster_size, num_highest_trans=self.num_highest_transparent_cluster, chi2=self.chi2)

		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.run1))
		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.run2))
		RecreateSoftLink(self.out_dir + '/' + self.subdir + '/' + str(self.run1), self.scratch_path, str(self.run1) + '_' + self.subdir, 'dir', False)
		RecreateSoftLink(self.out_dir + '/' + self.subdir + '/' + str(self.run2), self.scratch_path, str(self.run2) + '_' + self.subdir, 'dir', False)
		self.Print_subprocess_command('{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run1), self.settings_dir + '/' + self.subdir, self.out_dir + '/' + self.subdir, self.data_dir + '/' + str(self.run1))
		self.Print_subprocess_command('{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run2), self.settings_dir + '/' + self.subdir, self.out_dir + '/' + self.subdir, self.data_dir + '/' + str(self.run2))
		with open(os.devnull, 'w') as FNULL:
			self.sub_pro1 = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run1), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i',
		                                self.data_dir + '/' + str(self.run1)], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True)
			if self.batch:
				self.sub_pro2 = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run2), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i',
			                                self.data_dir + '/' + str(self.run2)], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True)
			else:
				self.sub_pro2 = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run2), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i',
		                                    self.data_dir + '/' + str(self.run2)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		while (self.sub_pro1.poll() is None) or (self.sub_pro2.poll() is None):
			time.sleep(5)
		if self.sub_pro1.poll() == 0 and self.sub_pro2.poll() == 0:
			print 'Run finished successfully'
		else:
			print 'Run could have failed. Obtained return code:', self.sub_pro1.poll(), 'and', self.sub_pro2.poll()
		CloseSubprocess(self.sub_pro1, True, False)
		CloseSubprocess(self.sub_pro2, True, False)

	def Merge_Files(self):
		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.merged_run_number))
		CreateDirectoryIfNecessary(self.data_dir + '/' + str(self.merged_run_number))
		if os.path.isfile('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number)):
			tempf = ro.TFile('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number))
			tempt = tempf.Get('rawTree')
			if tempt.GetEntries() >= self.total_events1 + self.total_events2:
				print 'File already merged :), continueing with the analysis :D...'
				self.run = self.merged_run_number
				self.total_events = self.total_events1 + self.total_events2
				self.num_events = self.num_events_merged
				self.num_events = self.total_events if self.num_events == 0 else self.num_events
				return
			else:
				tempf.Close()
		command = ['hadd', '{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number)]
		if os.path.isfile('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.run1)) and os.path.isfile('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.run2)):
			command.append('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.run1))
			command.append('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.run2))
		else:
			print 'At least one of the files corresponding for runs {r1} and {r2} does not exist! Exiting...'.format(d=self.out_dir, r1=self.run1, r2=self.run2)
			sys.exit(os.EX_SOFTWARE)
		print 'Merging files with the following command:', command, '...', ; sys.stdout.flush()
		self.sub_pro = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
		while self.sub_pro.poll() is None:
			time.sleep(5)
		print 'Done'
		CloseSubprocess(self.sub_pro, True, True)
		print 'Correcting event numbers...'
		tempf = ro.TFile('{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number), 'READ')
		tempt = tempf.Get('rawTree')
		event_new = np.zeros(1, 'I')
		tempt.SetBranchAddress('EventNumber', event_new)
		tempf_n = ro.TFile('{d}/{s}/{r}/rawData_n.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number), 'RECREATE')
		tempt_n = tempt.CloneTree(0)
		self.total_events = self.total_events1 + self.total_events2
		temp_bar = CreateProgressBarUtils(self.total_events)
		temp_bar.start()
		for ev in range(self.total_events):
			tempt.GetEntry(ev)
			event_new.fill(ev)
			tempt_n.Fill()
			event_new.fill(0)
			temp_bar.update(ev + 1)

		tempt_n.AutoSave()
		del tempf, tempf_n
		print 'blaaa'
		ipdb.set_trace()
		shutil.move('{d}/{s}/{r}/rawData_n.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number), '{d}/{s}/{r}/rawData.{r}.root'.format(d=self.out_dir, s=self.subdir, r=self.merged_run_number))
		self.run = self.merged_run_number
		self.num_events = self.num_events_merged
		self.num_events = self.total_events if self.num_events == 0 else self.num_events

	def Create_Run_List(self, do_single_ch=False):
		CreateDirectoryIfNecessary(self.run_lists_dir)
		ped = 1 if self.do_pedestal else 0
		clu = 1 if self.do_cluster else 0
		sele = 1 if self.do_selection else 0
		alig = 1 if self.do_alignment else 0
		tran = 1 if self.do_transparent else 0
		CreateDirectoryIfNecessary(self.run_lists_dir+'/{f}'.format(f=self.subdir))
		with open(self.run_lists_dir + '/{f}/RunList_'.format(f=self.subdir)+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t{al}\t0\t{t}\n#\n'.format(r=self.run, n=self.total_events, p=ped, c=clu, s=sele, al=alig, t=tran))
		if self.do_even or self.do_odd:
			CreateDirectoryIfNecessary(self.run_lists_dir+'/'+self.subdir+'/odd')
			with open(self.run_lists_dir+'/'+self.subdir+'/odd' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.total_events, p=ped, c=clu, s=sele))
			CreateDirectoryIfNecessary(self.run_lists_dir+'/'+self.subdir+'/even')
			with open(self.run_lists_dir+'/'+self.subdir+'/even' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.total_events, p=ped, c=clu, s=sele))
		if do_single_ch:
			CreateDirectoryIfNecessary(self.run_lists_dir+'/'+self.subdir+'/channels')
			with open(self.run_lists_dir + '/' + self.subdir + '/channels/RunList_' + str(self.run) + '.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.total_events, p=ped, c=clu, s=sele))

	def Check_settings_file(self):
		if not os.path.isdir(self.settings_dir + '/' + self.subdir):
			os.makedirs(self.settings_dir + '/' + self.subdir)
		if not os.path.isfile(self.settings_dir + '/' + self.subdir + '/settings.{r}.ini'.format(r=self.run)):
			CreateDefaultSettingsFile(self.settings_dir + '/' + self.subdir, self.run, self.total_events, dut_name=self.dut_name, dut_volt=self.dut_volt, ev_ini=self.first_event, num_evs_ana=self.num_events, dia_input=self.dia_input, dia_sat=self.dia_saturation, max_trans_clust=self.max_transparent_cluster_size, num_highest_trans=self.num_highest_transparent_cluster, chi2=self.chi2)
		else:
			with open(self.settings_dir + '/' + self.subdir + '/settings.{r}.ini'.format(r=self.run), 'r') as f0:
				with open(self.settings_dir + '/' + self.subdir + '/settings.{r}.ini.temp'.format(r=self.run), 'w') as ftemp:
					lines = f0.readlines()
					lines_params = [line.split('=')[0].replace(' ', '') for line in lines]
					if 'event_start' not in lines_params:
						ftemp.write('event_start = {v}\n'.format(v=self.first_event))
					if 'max_transparent_cluster_size' not in lines_params:
						ftemp.write('max_transparent_cluster_size = {d}\n'.format(d=self.max_transparent_cluster_size))
					if 'dia_saturation' not in lines_params:
						ftemp.write('dia_saturation = {d}\n'.format(d=self.dia_saturation))
					if 'eta_corr_limit' not in lines_params:
						ftemp.write('eta_corr_limit = {v}\n'.format(v=self.eta_corr_limit))
					if 'Iter_Size' not in lines_params:
						ftemp.write('Iter_Size = {v}\n'.format(v=self.ped_buffer))
					if 'dia_input' not in lines_params:
						ftemp.write('dia_input = {d}\n'.format(d=self.dia_input))
					if 'diamondName' not in lines_params:
						ftemp.write('diamondName = {dn}\n'.format(dn=self.dut_name))
					if 'voltage' not in lines_params:
						ftemp.write('voltage = {v}\n'.format(v=self.dut_volt))
					if 'num_highest_transparent_cluster' not in lines_params:
						ftemp.write('num_highest_transparent_cluster = {d}\n'.format(d=self.num_highest_transparent_cluster))
					if 'alignment_chi2' not in lines_params:
						ftemp.write('alignment_chi2 = {ch}\n'.format(ch=self.chi2))
					if 'alignment_training_track_number' not in lines_params:
						ftemp.write('alignment_training_track_number = {e}\n'.format(e=int(self.alignment_events)))
					if '3dShortAnalysis' not in lines_params:
						ftemp.write('3dShortAnalysis = {v};\n'.format(v=int(self.do_3d)))
					if '3dLongAnalysis' not in lines_params:
						ftemp.write('3dLongAnalysis = {v};\n'.format(v=int(self.do_3d)))
					if '3dTransparentAnalysis' not in lines_params:
						ftemp.write('3dTransparentAnalysis = {v};\n'.format(v=int(self.do_3d)))
					if 'transparentChi2' not in lines_params and self.trans_chi2:
						ftemp.write('transparentChi2 = {tc};\n'.format(tc=self.trans_chi2))
					if 'pulse_height_di_max' not in lines_params:
						ftemp.write('pulse_height_di_max = {ph};\n'.format(ph=self.ph_dia_max))
					if 'pulse_height_num_bins' not in lines_params:
						ftemp.write('pulse_height_num_bins = {n}\n'.format(n=self.ph_dia_bins))
					if 'yOffset3D' not in lines_params:
						ftemp.write('yOffset3D = {v}\n'.format(v=self.yOffset))
					if 'TransparentAlignment' not in lines_params:
						ftemp.write('TransparentAlignment = {v}\n'.format(v=int(self.trans_align)))
					for line in lines:
						if line.startswith('runNo'):
							ftemp.write('runNo = {r}\n'.format(r=self.run))
						elif line.startswith('diamondName'):
							ftemp.write('diamondName = {dn}\n'.format(dn=self.dut_name))
						elif line.startswith('voltage'):
							ftemp.write('voltage = {v}\n'.format(v=self.dut_volt))
						elif line.startswith('Event'):
							ftemp.write('Events = {e}\n'.format(e=self.total_events))
						elif line.startswith('event_start'):
							ftemp.write('event_start = {ei}\n'.format(ei=self.first_event))
						elif line.startswith('Events_an'):
							ftemp.write('Events_ana = {ea}\n'.format(ea=self.num_events))
						elif line.startswith('max_transpa'):
							ftemp.write('max_transparent_cluster_size = {d}\n'.format(d=self.max_transparent_cluster_size))
						elif line.startswith('dia_saturat'):
							ftemp.write('dia_saturation = {d}\n'.format(d=self.dia_saturation))
						elif line.startswith('Iter_Size'):
							ftemp.write('Iter_Size = {v}\n'.format(v=self.ped_buffer))
						elif line.startswith('dia_inp'):
							ftemp.write('dia_input = {d}\n'.format(d=self.dia_input))
						elif line.startswith('num_highest_trans'):
							ftemp.write('num_highest_transparent_cluster = {d}\n'.format(d=self.num_highest_transparent_cluster))
						elif line.startswith('alignment_chi2'):
							if self.chi2:
								ftemp.write('alignment_chi2 = {ch}\n'.format(ch=self.chi2))
							else:
								ftemp.write(line)
						elif line.startswith('alignment_training_track_numb'):
							ftemp.write('alignment_training_track_number = {e}\n'.format(e=int(self.alignment_events)))
						elif line.startswith('3dShortAna'):
							ftemp.write('3dShortAnalysis = {v};\n'.format(v=int(self.do_3d)))
						elif line.startswith('3dLongAna'):
							ftemp.write('3dLongAnalysis = {v};\n'.format(v=int(self.do_3d)))
						elif line.startswith('3dTransparentAna'):
							ftemp.write('3dTransparentAnalysis = {v};\n'.format(v=int(self.do_3d)))
						elif line.startswith('TransparentAlignm'):
							ftemp.write('TransparentAlignment = {v}\n'.format(v=int(self.trans_align)))
						elif line.startswith('transparentChi2'):
							if self.trans_chi2:
								ftemp.write('transparentChi2 = {tc}\n'.format(tc=self.trans_chi2))
							else:
								ftemp.write(line)
						elif line.startswith('pulse_height_di_'):
							ftemp.write('pulse_height_di_max = {ph}\n'.format(ph=self.ph_dia_max))
						elif line.startswith('pulse_height_num_b'):
							ftemp.write('pulse_height_num_bins = {nb}\n'.format(nb=self.ph_dia_bins))
						elif line.startswith('yOffset3'):
							if self.yOffset != 0:
								ftemp.write('yOffset3D = {v}\n'.format(v=self.yOffset))
							else:
								ftemp.write(line)
						elif line.startswith('eta_corr_lim'):
							ftemp.write('eta_corr_limit = {v}\n'.format(v=self.eta_corr_limit))
						else:
							ftemp.write(line)
			shutil.move(self.settings_dir + '/' + self.subdir + '/settings.{r}.ini.temp'.format(r=self.run), self.settings_dir + '/' + self.subdir + '/settings.{r}.ini'.format(r=self.run))

	def CheckStripTelescopeAnalysis(self):
		if os.path.isdir(self.StripTelescopeAnalysis_path):
			if not os.path.isfile(self.StripTelescopeAnalysis_path + '/diamondAnalysis'):
				ExitMessage('{p}/diamondAnalysis does not exist. Exiting'.format(p=self.StripTelescopeAnalysis_path))
		else:
			ExitMessage('{d} does not exist. Exiting'.format(d=self.StripTelescopeAnalysis_path))

	def First_Analysis(self):
		self.subdir = 'no_mask'
		if self.delete_old:
			self.Delete_old()
		print 'Starting first analysis (no_mask)...'
		self.RunAnalysis()
		print 'Finished with first analysis :)'

	def RunAnalysis(self):
		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.run))
		RecreateSoftLink(self.out_dir + '/' + self.subdir + '/' + str(self.run), self.scratch_path, str(self.run) + '_' + self.subdir, 'dir', False)
		self.Print_subprocess_command('{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), self.settings_dir + '/' + self.subdir, self.out_dir + '/' + self.subdir, self.data_dir + '/' + str(self.run))
		if self.batch:
			print 'No output'
		else:
			print 'With verbose'
		with open(os.devnull, 'w') as FNULL:
			command = ['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir, '-o', self.out_dir + '/' + self.subdir, '-i', self.data_dir + '/' + str(self.run)]
			if self.do_3d:
				command.append('-3d')
				command.append('1')
			print command
			if self.batch:
				self.sub_pro = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, stdout=FNULL, stderr=subp.STDOUT,close_fds=True)
			else:
				self.sub_pro = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, close_fds=True)
			while self.sub_pro.poll() is None:
				time.sleep(2)
			if self.sub_pro.poll() == 0:
				print 'Run', self.run, 'finished successfully'
			else:
				print 'Run', self.run, 'could have failed. Obtained return code:', self.sub_pro.poll()
			CloseSubprocess(self.sub_pro, True, False)

	def Print_subprocess_command(self, runlist, setting, outdir, inputdir):
		print 'Executing:\n{p}/diamondAnalysis -r {r} -s {s} -o {o} -i {i}\n'.format(p=self.StripTelescopeAnalysis_path, r=runlist, s=setting, o=outdir, i=inputdir)

	def LinkRootFiles(self, source, dest, upto='selection', doSymlink=True, doCopy=True, nodir=False):
		steps = ['raw', 'pedestal', 'cluster', 'selection', 'align', 'transparent', '3d']
		cumulative = []
		for elem in steps:
			cumulative.append(elem)
			if elem == upto:
				break
		successful = True
		if 'raw' in cumulative:
			if not os.path.isfile(dest + '/rawData.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/rawData.{r}.root'.format(r=self.run), dest, 'rawData.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
		if 'pedestal' in cumulative:
			if not os.path.isfile(dest + '/pedestalData.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/pedestalData.{r}.root'.format(r=self.run), dest, 'pedestalData.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
				if not nodir: RecreateLink(source + '/pedestalAnalysis', dest, 'pedestalAnalysis', doSymlink, doCopy)
		if 'cluster' in cumulative:
			if not os.path.isfile(dest + '/clusterData.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/clusterData.{r}.root'.format(r=self.run), dest, 'clusterData.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
				if not nodir: RecreateLink(source + '/clustering', dest, 'clustering', doSymlink, doCopy)
			if not os.path.isfile(dest + '/etaCorrection.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/etaCorrection.{r}.root'.format(r=self.run), dest, 'etaCorrection.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
			if os.path.isfile(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run)) and doCopy:
				shutil.copy(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run), dest + '/')
		if 'selection' in cumulative:
			if not os.path.isfile(dest + '/selectionData.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/selectionData.{r}.root'.format(r=self.run), dest, 'selectionData.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
				if not nodir: RecreateLink(source + '/selectionAnalysis', dest, 'selectionAnalysis', doSymlink, doCopy)
				if not nodir: RecreateLink(source + '/selections', dest, 'selections', doSymlink, doCopy)
		if 'align' in cumulative:
			if not os.path.isfile(dest + '/alignment.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/alignment.{r}.root'.format(r=self.run), dest, 'alignment.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
				if not nodir: RecreateLink(source + '/alignment', dest, 'alignment', doSymlink, doCopy)
		if 'transparent' in cumulative:
			if not os.path.isfile(dest + '/transparentAnalysis'):
				if not nodir:
					successful2 = RecreateLink(source + '/transparentAnalysis', dest, 'transparentAnalysis', doSymlink, doCopy)
					successful = successful and successful2
		if '3d' in cumulative:
			if not os.path.isfile(dest + '/analysis3d.root.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/analysis3d.root.{r}.root'.format(r=self.run), dest, 'analysis3d.root.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
			if not os.path.isfile(dest + '/analysis3d-2.root.{r}.root'.format(r=self.run)):
				successful2 = RecreateLink(source + '/analysis3d-2.root.{r}.root'.format(r=self.run), dest, 'analysis3d-2.root.{r}.root'.format(r=self.run), doSymlink, doCopy)
				successful = successful and successful2
			if not os.path.isfile(dest + '/3dDiamondAnalysis'):
				if not nodir:
					successful2 = RecreateLink(source + '/3dDiamondAnalysis', dest, '3dDiamondAnalysis', doSymlink, doCopy)
					successful = successful and successful2

		if os.path.isfile(source + '/index.php'):
			shutil.copyfile(source + '/index.php', dest + '/index.php')
		if os.path.isfile(source + '/overview.html'):
			shutil.copyfile(source + '/overview.html', dest + '/overview.html')
		if os.path.isfile(source + '/results_{r}.res'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.res'.format(r=self.run), dest + '/results_{r}.res'.format(r=self.run))
		if os.path.isfile(source + '/results_{r}.txt'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.txt'.format(r=self.run), dest + '/results_{r}.txt'.format(r=self.run))
		if not os.path.isfile(dest + '/Results.{r}.root'.format(r=self.run)):
			RecreateLink(source + '/Results.{r}.root'.format(r=self.run), dest, 'Results.{r}.root'.format(r=self.run), doSymlink, doCopy)
		return successful

	def Normal_Analysis(self):
		if self.delete_old:
			self.Delete_old()
		if os.path.isfile(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
			if self.Get_num_events(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'rawTree') < self.num_events:
				print 'Rerunning first event (no_mask) as it has less entries than required...'
				self.subdir = 'no_mask'
				self.num_events = self.total_events
				self.Create_Run_List()
				self.Check_settings_file()
				self.First_Analysis()
			elif self.Get_num_events(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'rawTree') > self.num_events:
				print 'Need to recreate the rawTree'
				self.ExtractFromOriginalRawTree(self.subdir)
		# elif os.path.isfile(self.out_dir + '/no_mask/{r}/rawData.{r}.root'.format(r=self.run)):
		# 	if self.Get_num_events(self.out_dir + '/no_mask/{r}/rawData.{r}.root'.format(r=self.run), 'rawTree') <= self.num_events:
		# 		# print 'Extracting', self.num_events, 'events from no_mask run'
		# 		# self.ExtractFromOriginalRawTree('no_mask')
		# 	# else:
		# 		if self.Get_num_events(self.out_dir + '/no_mask/{r}/rawData.{r}.root'.format(r=self.run), 'rawTree') < self.num_events:
		# 			print 'The no_mask raw file has', self.Get_num_events(self.out_dir + '/no_mask/{r}/rawData.{r}.root'.format(r=self.run), 'rawTree'), 'events which is less than the requested. Will analyse all events'
		# 			self.num_events = self.Get_num_events(self.out_dir + '/no_mask/{r}/rawData.{r}.root'.format(r=self.run), 'rawTree')
		# 		self.LinkRootFiles(self.out_dir + '/no_mask', self.out_dir + '/' + self.subdir, 'raw', True, True)
		self.RunAnalysis()
		if self.do_cross_talk:
			if os.path.isfile(self.cross_talk_correction_path):
				self.current_dir = self.out_dir + '/' + self.subdir + '/' + str(self.run)
				os.chdir(self.current_dir)
				lines = []
				sil_list = []
				sil_value = 0
				dia = 0
				if os.path.isfile('crossTalkCorrectionFactors.{r}.txt'.format(r=self.run)):
					with open('crossTalkCorrectionFactors.{r}.txt'.format(r=self.run), 'r') as f0:
						lines = f0.readlines()
					for line in lines:
						if int(line[0])<8:
							sil_list.append(float(line.split(': ')[1].split('%')[0]))
						elif int(line[0]) == 8:
							dia = float(line.split(': ')[1].split('%')[0])
					sil_values = np.array(sil_list, 'f8')
					sil_mean = sil_values.mean(dtype='f8')

					print 'Running from run outpu directory:', self.cross_talk_correction_path, '-o', '.', '-r', str(self.run), '-s', str(sil_mean), '-s0', str(sil_values[0]), '-s1', str(sil_values[1]), '-s2', str(sil_values[2]), '-s3', str(sil_values[3]), '-s4', str(sil_values[4]), '-s5', str(sil_values[5]), '-s6', str(sil_values[6]), '-s7', str(sil_values[7]), '-d', str(dia), '-ss', str(255), '-ds', str(self.dia_saturation)

					with open(os.devnull, 'w') as FNULL:
						if self.batch:
							self.sub_pro_crosstalk = subp.Popen([self.cross_talk_correction_path, '-o', '.', '-r', str(self.run), '-s', str(sil_mean), '-s0', str(sil_values[0]), '-s1', str(sil_values[1]), '-s2', str(sil_values[2]), '-s3', str(sil_values[3]), '-s4', str(sil_values[4]), '-s5', str(sil_values[5]), '-s6', str(sil_values[6]), '-s7', str(sil_values[7]), '-d', str(dia), '-ss', str(255), '-ds', str(self.dia_saturation)], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, stderr=subp.STDOUT, close_fds=True)
						else:
							self.sub_pro_crosstalk = subp.Popen([self.cross_talk_correction_path, '-o', '.', '-r', str(self.run), '-s', str(sil_mean), '-s0', str(sil_values[0]), '-s1', str(sil_values[1]), '-s2', str(sil_values[2]), '-s3', str(sil_values[3]), '-s4', str(sil_values[4]), '-s5', str(sil_values[5]), '-s6', str(sil_values[6]), '-s7', str(sil_values[7]), '-d', str(dia), '-ss', str(255), '-ds', str(self.dia_saturation)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
						while self.sub_pro_crosstalk.poll() is None:
							time.sleep(3)
						if self.sub_pro_crosstalk.poll() == 0:
							print 'Finished correcting the raw faile successfully :)'
						else:
							print 'The correction could have failed. Obtained return code:', self.sub_pro_crosstalk.poll()
						CloseSubprocess(self.sub_pro_crosstalk, True, False)
					# ro.gROOT.ProcessLine('.x {p}/createAsymmetricEtaSample.C({r},{s},{d},-1)'.format(p=self.StripTelescopeAnalysis_path, r=self.run, s=sil_mean, d=dia))
					if not os.path.isdir('../../cross_{d}/{r}'.format(d=self.subdir, r=self.run)):
						os.makedirs('../../cross_{d}/{r}'.format(d=self.subdir, r=self.run))
					self.current_dir = self.out_dir + '/cross_' + self.subdir + '/' + str(self.run)
					os.chdir(self.current_dir)
					if self.symlinks:
						try:
							os.symlink('../../{d}/{r}/rawData.{r}0.root'.format(d=self.subdir, r=self.run), 'rawData.{r}.root'.format(r=self.run))
						except Exception:
							print 'Could not ln -s the corrected raw file for further anlysis. Do it manually.'
					else:
						try:
							shutil.copy('../../{d}/{r}/{n}'.format(d=self.subdir, r=self.run, n=os.readlink('../../{d}/{r}/rawData.{r}0.root'.format(d=self.subdir, r=self.run))), 'rawData.{r}.root'.format(r=self.run))
						except Exception:
							print 'Could not copy the corrected raw file for further anlysis. Do it manually.'
					self.current_dir = self.working_dir
					os.chdir(self.working_dir)
					if not os.path.isdir('{d}/cross_{sd}'.format(d=self.settings_dir, sd=self.subdir)):
						os.makedirs('{d}/cross_{sd}'.format(d=self.settings_dir, sd=self.subdir))
					try:
						shutil.copy('{d}/{sd}/settings.{r}.ini'.format(d=self.settings_dir, sd=self.subdir, r=self.run), '{d}/cross_{sd}/settings.{r}.ini'.format(d=self.settings_dir, sd=self.subdir, r=self.run))
						print 'Copied settings file to: {d}/{sd}/settings.{r}.ini'.format(d=self.settings_dir, sd=self.subdir, r=self.run), '. Edit it if you want for the analysis with the feed-across corrected raw file.'
					except Exception:
						print 'Could not copy the settings file for further anlysis of the feed-across corrected data. Do it manually.'

				else:
					print "The crossTalkCorrectionFactors text file does not exist. Can't run the correction"
			else:
				print "the path", self.cross_talk_correction_path, "does not exist. Can't run the correction"
		print 'Finished :)'

	def Delete_old(self, upto='3d'):
		print 'Deleting old upto', upto, ; sys.stdout.flush()
		steps = ['raw', 'pedestal', 'cluster', 'selection', 'align', 'transparent', '3d']
		cumulative = []
		for elem in steps:
			cumulative.append(elem)
			if elem == upto:
				break
		stem_dir = '{od}/{sd}/{r}'.format(od=self.out_dir, sd=self.subdir, r=self.run)
		if 'raw' in cumulative:
			if os.path.isfile('{sd}/rawData.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/rawData.{r}.root'.format(sd=stem_dir, r=self.run))
		if 'pedestal' in cumulative:
			if os.path.isfile('{sd}/pedestalData.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/pedestalData.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.islink('{sd}/pedestalAnalysis'.format(sd=stem_dir)):
				os.unlink('{sd}/pedestalAnalysis'.format(sd=stem_dir))
			elif os.path.isdir('{sd}/pedestalAnalysis'.format(sd=stem_dir)):
				shutil.rmtree('{sd}/pedestalAnalysis'.format(sd=stem_dir))
		if 'cluster' in cumulative:
			if os.path.isfile('{sd}/clusterData.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/clusterData.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.isfile('{sd}/etaCorrection.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/etaCorrection.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.isfile('{sd}/crossTalkCorrectionFactors.{r}.txt'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/crossTalkCorrectionFactors.{r}.txt'.format(sd=stem_dir, r=self.run))
			if os.path.islink('{sd}/clustering'.format(sd=stem_dir)):
				os.unlink('{sd}/clustering'.format(sd=stem_dir))
			elif os.path.isdir('{sd}/clustering'.format(sd=stem_dir)):
				shutil.rmtree('{sd}/clustering'.format(sd=stem_dir))
		if 'align' in cumulative:
			if os.path.isfile('{sd}/alignment.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/alignment.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.islink('{sd}/alignment'.format(sd=stem_dir)):
				os.unlink('{sd}/alignment'.format(sd=stem_dir))
			elif os.path.isdir('{sd}/alignment'.format(sd=stem_dir)):
				shutil.rmtree('{sd}/alignment'.format(sd=stem_dir))
		if 'transparent' in cumulative:
			if os.path.islink('{sd}/transparentAnalysis'.format(sd=stem_dir)):
				os.unlink('{sd}/transparentAnalysis'.format(sd=stem_dir))
			elif os.path.isdir('{sd}/transparentAnalysis'.format(sd=stem_dir)):
				shutil.rmtree('{sd}/transparentAnalysis'.format(sd=stem_dir))
		if '3d' in cumulative:
			if os.path.isfile('{sd}/analysis3d.root.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/analysis3d.root.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.isfile('{sd}/analysis3d-2.root.{r}.root'.format(sd=stem_dir, r=self.run)):
				os.unlink('{sd}/analysis3d-2.root.{r}.root'.format(sd=stem_dir, r=self.run))
			if os.path.islink('{sd}/3dDiamondAnalysis'.format(sd=stem_dir)):
				os.unlink('{sd}/3dDiamondAnalysis'.format(sd=stem_dir))
			elif os.path.isdir('{sd}/3dDiamondAnalysis'.format(sd=stem_dir)):
				shutil.rmtree('{sd}/3dDiamondAnalysis'.format(sd=stem_dir))
		print 'Done'

	def Get_num_events(self, rootfile, treename):
		tempf = ro.TFile(rootfile, 'READ')
		tempt = tempf.Get(treename)
		num_evts = tempt.GetEntries()
		tempf.Close()
		return num_evts

	def ExtractFromOriginalRawTree(self, originalsubdir='no_mask'):
		CreateDirectoryIfNecessary(self.out_dir + '/' + originalsubdir + '/{r}'.format(r=self.run))
		if os.path.isfile(self.out_dir + '/' + originalsubdir + '/{r}/rawData.{r}.root'.format(r=self.run)):
			tempf = ro.TFile(self.out_dir + '/' + originalsubdir + '/{r}/rawData.{r}.root'.format(r=self.run), 'READ')
			tempt = tempf.Get('rawTree')
			print 'Extracting only', self.num_events, 'events starting from', self.first_event, 'to analyse...', ; sys.stdout.flush()
			leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_event, eva=self.num_events))
			while leng > tempt.GetEstimate():
				tempt.SetEstimate(leng)
				leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_event, eva=self.num_events))
			evlist = ro.gDirectory.Get('evlist')
			tempt.SetEventList(evlist)
			tempnf = ro.TFile(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData2.{r}.root'.format(r=self.run), 'RECREATE')
			tempnt = tempt.CopyTree('')
			tempnt.SetName('rawTree')
			tempnt.Write()
			tempnf.Close()
			tempf.Close()
			if os.path.islink(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
				os.unlink(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run))
			shutil.move(self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData2.{r}.root'.format(r=self.run), self.out_dir + '/' + self.subdir + '/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run))
			print 'Done'
		else:
			ExitMessage('Cannot extract from ' + originalsubdir + ' as it does not exist. Run first the analysis with sub directory: ' + originalsubdir)

	def GetIndividualChannelHitmap(self):
		print 'Starting channel occupancy plots...'
		cont = False
		CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep'.format(sd=self.subdir, r=self.run))
		for ch in xrange(diaChs):
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=self.subdir, c=ch, r=self.run))
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=self.subdir, c=ch, r=self.run))
			do_continue = self.LinkRootFiles(self.out_dir + '/'+self.subdir+'/' + str(self.run), self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=self.subdir, c=ch, r=self.run), 'cluster', doSymlink=True, doCopy=False, nodir=True)
			if not do_continue:
				ExitMessage('Cannot create symlinks. Exiting...')
			self.ModifySettingsOneChannel(ch, self.subdir)
		num_cores = mp.cpu_count()
		num_parrallel = 1 if num_cores < 4 else 2 if num_cores < 8 else 4 if num_cores < 16 else 8 if num_cores < 32 else 16
		num_jobs = diaChs / num_parrallel
		channel_runs = np.arange(diaChs).reshape([num_jobs, num_parrallel])
		for bat in xrange(channel_runs.shape[0]):
			procx = []
			for proc in xrange(channel_runs.shape[1]):
				print 'Getting channel:', channel_runs[bat, proc], '...'
				procx.append(subp.Popen(
					['diamondAnalysis', '-r', '{d}/{sd}/channels/RunList_{r}.ini'.format(sd=self.subdir, d=self.run_lists_dir, r=self.run), '-s', self.settings_dir + '/{sd}/channels/{c}'.format(sd=self.subdir, c=channel_runs[bat, proc]),
					 '-o', self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=self.subdir, c=channel_runs[bat, proc], r=self.run), '-i', self.data_dir], bufsize=-1, stdin=subp.PIPE,
					stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True))
			for proc in xrange(channel_runs.shape[1]):
				while procx[proc].poll() is None:
					time.sleep(3)
				CloseSubprocess(procx[proc], stdin=True, stdout=False)
				print 'Done with channel:', channel_runs[bat, proc], ':)'
			del procx

	def ModifySettingsOneChannel(self, ch=0, sub_dir='no_mask'):
		CreateDirectoryIfNecessary(self.settings_dir + '/{sd}/channels'.format(sd=sub_dir))
		CreateDirectoryIfNecessary(self.settings_dir + '/{sd}/channels/{c}'.format(sd=sub_dir, c=ch))
		with open(self.settings_dir + '/{sd}/settings.{r}.ini'.format(sd=sub_dir, r=self.run), 'r') as fin:
			with open(self.settings_dir + '/{sd}/channels/{c}/settings.{r}.ini'.format(sd=sub_dir, c=ch, r=self.run), 'w') as fch:
				for line in fin:
					if not line.startswith('Dia_channel_screen_channels'):
						fch.write(line)
					else:
						channel_str_new = self.ModifyStringOneChannel(ch)
						fch.write('Dia_channel_screen_channels = {' + channel_str_new + '}\n')

	def ModifyStringOneChannel(self, ch=0):
		if ch == 0:
			return '1-{d}'.format(d=diaChs - 1)
		elif ch == diaChs - 1:
			return '0-{d}'.format(d=diaChs - 2)
		else:
			return '0-{cp},{cn}-{d}'.format(cp=ch - 1, cn = ch + 1, d=diaChs - 1)

def main():
	parser = OptionParser()
	parser.add_option('-w', '--workingdir', dest='workingdir', type='string', help='Working directory')
	parser.add_option('-s', '--settings', dest='settings_f', type='string', help='Settings file containing the information on the run and the analysis to do (e.g. settings.ini)')
	parser.add_option('--first', dest='first', default=False, action='store_true', help='enables first analysis wich has everything un-masked')
	parser.add_option('--normal', dest='normal', default=False, action='store_true', help='enables normal analysis')
	parser.add_option('--singlechannel', dest='singlech', default=False, action='store_true', help='enables single channel study. Requires a preexisting first analysis')
	parser.add_option('-q', '--quiet', dest='quiet', default=False, action='store_true', help='enables quiet mode: no verbose')

	(options, args) = parser.parse_args()
	working_dir = str(options.workingdir)
	settings_f = str(options.settings_f)
	first_ana = bool(options.first)
	normal_ana = bool(options.normal)
	single_ch = bool(options.singlech)
	verb = not bool(options.quiet)

	rd42 = RD42Analysis(working_dir=working_dir, verb=verb)
	rd42.ReadInputFile(settings_f)
	if rd42.do_merged:
		rd42.Convert_Files()
		rd42.Merge_Files()
	if first_ana:
		rd42.subdir = 'no_mask'
	rd42.Create_Run_List(do_single_ch=single_ch)
	rd42.Check_settings_file()
	rd42.CheckStripTelescopeAnalysis()

	if first_ana:
		print 'Starting first analysis (no_mask)...\n'
		rd42.First_Analysis()
	elif normal_ana:
		print 'Starting normal analysis...\n'
		rd42.Normal_Analysis()
	if single_ch:
		rd42.GetIndividualChannelHitmap()


if __name__ == '__main__':
	main()
