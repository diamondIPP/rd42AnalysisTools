#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
from ConfigParser import ConfigParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
from NoiseExtraction import NoiseExtraction
import os, sys, shutil
from Utils import *
import subprocess as subp
import multiprocessing as mp

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
		self.symlinks = True
		# self.scratch_path = '/scratch/strip_telescope_tests/runDiego/output'  # at snickers

		self.delete_old = False
		self.first_event = 0
		self.num_events = 0
		self.do_pedestal = False
		self.do_cluster = False
		self.do_selection = False
		self.do_alignment = False
		self.do_transparent = False
		self.do_3d = False

		self.sub_pro, self.sub_pro_e, self.sub_pro_o = None, None, None
		self.process_f = None
		self.process_e = None
		self.process_o = None
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
					if pars.has_option('RUN', 'symlinks'):
						self.symlinks = pars.getboolean('RUN', 'symlinks')
					if pars.has_option('RUN', 'delete_old'):
						self.delete_old = pars.getboolean('RUN', 'delete_old')

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
		self.Modify_even_odd()

	def Modify_even_odd(self):
		self.Modify_even()
		self.Modify_odd()

	def Modify_even(self):
		print 'Modifying even settings file...', ; sys.stdout.flush()
		Replace_Settings_Line(self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run), 'Dia_channel_screen_channels', 'even')
		print 'Done'

	def Modify_odd(self):
		print 'Modifying odd settings file...', ; sys.stdout.flush()
		Replace_Settings_Line(self.settings_dir + '/odd/settings.{r}.ini'.format(r=self.run), 'Dia_channel_screen_channels', 'odd')
		print 'Done'

	def CheckStripTelescopeAnalysis(self):
		if os.path.isdir(self.StripTelescopeAnalysis_path):
			if not os.path.isfile(self.StripTelescopeAnalysis_path + '/diamondAnalysis'):
				ExitMessage('{p}/diamondAnalysis does not exist. Exiting'.format(p=self.StripTelescopeAnalysis_path))
		else:
			ExitMessage('{d} does not exist. Exiting'.format(d=self.StripTelescopeAnalysis_path))

	def Print_subprocess_command(self, runlist, setting, outdir, inputdir):
		print 'Executing:\n{p}/diamondAnalysis -r {r} -s {s} -o {o} -i {i}\n'.format(p=self.StripTelescopeAnalysis_path, r=runlist, s=setting, o=outdir, i=inputdir)

	def RunAnalysis(self):
		CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/' + str(self.run))
		if self.symlinks:
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
			self.sub_pro_o = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/odd/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir + '/odd', '-o', self.out_dir + '/' + self.subdir + '/odd', '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), close_fds=True)
		if self.do_even:
			CreateDirectoryIfNecessary(self.out_dir + '/' + self.subdir + '/even/' + str(self.run))
			self.Print_subprocess_command('{d}/{sd}/even/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), self.settings_dir + '/' + self.subdir + '/even', self.out_dir + '/' + self.subdir + '/even', self.data_dir + '/' + str(self.run))
			self.sub_pro_e = subp.Popen(['{p}/diamondAnalysis'.format(p=self.StripTelescopeAnalysis_path), '-r', '{d}/{sd}/even/RunList_{r}.ini'.format(d=self.run_lists_dir, sd=self.subdir, r=self.run), '-s', self.settings_dir + '/' + self.subdir + '/even', '-o', self.out_dir + '/' + self.subdir + '/even', '-i', self.data_dir + '/' + str(self.run)], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), close_fds=True)
		if self.do_odd:
			while self.sub_pro_o.poll() is None:
				pass
			if self.sub_pro_o.poll() == 0:
				print 'Run odd finished'
			else:
				print 'Run odd could have failed. Obtained return code:', self.sub_pro_o.poll()
			CloseSubprocess(self.sub_pro_o, True, False)
		if self.do_even:
			while self.sub_pro_e.poll() is None:
				pass
			if self.sub_pro_e.poll() == 0:
				print 'Run even finished'
			else:
				print 'Run even could have failed. Obtained return code:', self.sub_pro_e.poll()
			CloseSubprocess(self.sub_pro_e, True, False)

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

	def First_Analysis(self):
		self.subdir = 'no_mask'
		if self.delete_old:
			self.Delete_old()
		print 'Starting first analysis (no_mask)...'
		self.RunAnalysis()
		print 'Finished with first analysis :)'

	def GetIndividualChannelHitmap(self):
		print 'Starting channel occupancy plots...'
		cont = False
		CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep'.format(sd=self.subdir, r=self.run))
		for ch in xrange(diaChs):
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=self.subdir, c=ch, r=self.run))
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=self.subdir, c=ch, r=self.run))
			if self.symlinks:
				self.LinkRootFiles(self.out_dir + '/{sd}/' + str(self.run), self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=self.subdir, c=ch, r=self.run), 'cluster')
			else:
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
					['diamondAnalysis', '-r', '{d}/no_mask/RunList_{r}.ini'.format(sd=self.subdir, d=self.run_lists_dir, r=self.run), '-s', self.settings_dir + '/{sd}/channels/{c}'.format(sd=self.subdir, c=channel_runs[bat, proc]),
					 '-o', self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=self.subdir, c=channel_runs[bat, proc], r=self.run), '-i', self.data_dir], bufsize=-1, stdin=subp.PIPE,
					stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True))
			for proc in xrange(channel_runs.shape[1]):
				while procx[proc].poll() is None:
					# procx[proc].stdout.
					continue
			for proc in xrange(channel_runs.shape[1]):
				CloseSubprocess(procx[proc], stdin=True, stdout=False)
				print 'Done with channel:', channel_runs[bat, proc], ':)'
			del procx

	def Normal_Analysis(self):
		if self.delete_old:
			self.Delete_old()
		if os.path.isfile(self.out_dir + '/no_mask/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
			if self.first_event == 0 and self.num_events == self.total_events:
				print 'Linking raw file from no_mask...',; sys.stdout.flush()
				self.LinkRootFiles(self.out_dir + '/no_mask/' + str(self.run), self.out_dir + '/full/' + str(self.run), 'raw')
				print 'Done'
			else:
				recreate = True
				if os.path.isfile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)) or os.path.islink(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
					tempf0 = ro.TFile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'READ')
					tempt0 = tempf0.Get('rawTree')
					recreate = False if tempt0.GetEntries() == self.num_events else True
					tempf0.Close()
				if recreate:
					tempf = ro.TFile(self.out_dir + '/no_mask/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'READ')
					tempt = tempf.Get('rawTree')
					print 'Extracting only {eva} events starting from {evi} for analysis...'.format(eva=self.num_events, evi=self.first_event),; sys.stdout.flush()
					leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_event, eva=self.num_events))
					# leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{evf}-{evi})<=({evf}-{evi})'.format(evi=self.first_event, evf=self.first_event+self.num_events-1))
					while leng > tempt.GetEstimate():
						tempt.SetEstimate(leng)
						leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_event, eva=self.num_events))
					# leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{evf}-{evi})<=({evf}-{evi})'.format(evi=self.first_event, evf=self.first_event+self.num_events-1))
					evlist = ro.gDirectory.Get('evlist')
					tempt.SetEventList(evlist)
					tempnf = ro.TFile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'RECREATE')
					tempnt = tempt.CopyTree('')
					tempnt.Write()
					tempnf.Close()
					tempf.Close()
					print 'Done'
		self.process_f = subp.Popen(['diamondAnalysis', '-r', '{d}/full/RunList_{r}.ini'.format(d=self.run_lists_dir, r=self.run), '-s', self.settings_dir + '/full', '-o', self.out_dir + '/full', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		is_finished_cluster = False
		is_even_odd_started = False
		while self.process_f.poll() is None:
			if not is_even_odd_started:
				is_finished_cluster = self.Check_if_clustering_finished()
				if is_finished_cluster:
					self.LinkRootFiles(self.out_dir + '/full/' + str(self.run), self.out_dir + '/even/' + str(self.run), 'cluster', doCopy=False)
					self.LinkRootFiles(self.out_dir + '/full/' + str(self.run), self.out_dir + '/odd/' + str(self.run), 'cluster', doCopy=False)
					self.process_e = subp.Popen(['diamondAnalysis', '-r', '{d}/even/RunList_{r}.ini'.format(d=self.run_lists_dir, r=self.run), '-s', self.settings_dir + '/even', '-o', self.out_dir + '/even', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True)
					self.process_o = subp.Popen(['diamondAnalysis', '-r', '{d}/odd/RunList_{r}.ini'.format(d=self.run_lists_dir, r=self.run), '-s', self.settings_dir + '/odd', '-o', self.out_dir + '/odd', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True)
					is_even_odd_started = True
		print 'Full process finished'
		CloseSubprocess(self.process_f, True, False)
		temp_flag = True
		while self.process_e.poll() is None:
			if temp_flag:
				print 'waiting for even events to finish...'
				temp_flag = False
		CloseSubprocess(self.process_e, True, False)
		temp_flag = True
		while self.process_o.poll() is None:
			if temp_flag:
				print 'waiting for even events to finish...'
				temp_flag = False
		CloseSubprocess(self.process_o, True, False)
		print 'Finished :)'

	def Check_if_clustering_finished(self):
		if os.path.isdir(self.out_dir + '/full/{r}/selections'.format(r=self.run)):
			return True
		return False

	def LinkRootFiles(self, source, dest, upto='selection', doCopy=True):
		steps = ['raw', 'pedestal', 'cluster', 'selection', 'align', 'transparent']
		cummulative = []
		for elem in steps:
			cummulative.append(elem)
			if elem == upto:
				break
		if 'raw' in cummulative:
			if not os.path.isfile(dest + '/rawData.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/rawData.{r}.root'.format(r=self.run), dest, 'rawData.{r}.root'.format(r=self.run))
		if 'pedestal' in cummulative:
			if not os.path.isfile(dest + '/pedestalData.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/pedestalData.{r}.root'.format(r=self.run), dest, 'pedestalData.{r}.root'.format(r=self.run))
				RecreateLink(source + '/pedestalAnalysis', dest, 'pedestalAnalysis')
		if 'cluster' in cummulative:
			if not os.path.isfile(dest + '/clusterData.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/clusterData.{r}.root'.format(r=self.run), dest, 'clusterData.{r}.root'.format(r=self.run))
				RecreateLink(source + '/clustering', dest, 'clustering')
			if not os.path.isfile(dest + '/etaCorrection.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/etaCorrection.{r}.root'.format(r=self.run), dest, 'etaCorrection.{r}.root'.format(r=self.run))
			if os.path.isfile(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run)) and doCopy:
				shutil.copy(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run), dest + '/')
		if 'selection' in cummulative:
			if not os.path.isfile(dest + '/selectionData.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/selectionData.{r}.root'.format(r=self.run), dest, 'selectionData.{r}.root'.format(r=self.run))
				RecreateLink(source + '/selectionAnalysis', dest, 'selectionAnalysis')
				RecreateLink(source + '/selections', dest, 'selections')
		if 'align' in cummulative:
			if not os.path.isfile(dest + '/alignment.{r}.root'.format(r=self.run)):
				RecreateLink(source + '/alignment.{r}.root'.format(r=self.run), dest, 'alignment.{r}.root'.format(r=self.run))
				RecreateLink(source + '/alignment', dest, 'alignment')
		if 'transparent' in cummulative:
			if not os.path.isfile(dest + '/transparentAnalysis'):
				RecreateLink(source + '/transparentAnalysis', dest, 'transparentAnalysis')
		if os.path.isfile(source + '/index.php'):
			shutil.copyfile(source + '/index.php', dest + '/index.php')
		if os.path.isfile(source + '/overview.html'):
			shutil.copyfile(source + '/overview.html', dest + '/overview.html')
		if os.path.isfile(source + '/results_{r}.res'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.res'.format(r=self.run), dest + '/results_{r}.res'.format(r=self.run))
		if os.path.isfile(source + '/results_{r}.txt'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.txt'.format(r=self.run), dest + '/results_{r}.txt'.format(r=self.run))
		if not os.path.isfile(dest + '/Results.{r}.root'.format(r=self.run)):
			RecreateLink(source + '/Results.{r}.root'.format(r=self.run), dest, 'Results.{r}.root'.format(r=self.run))

	def RunNormalAnalysis(self, subdir=''):
		cont = False
		sub_dir = subdir
		valid_subdir_options = ['no_mask', 'full', '3D', 'strip', 'phantom', 'poly']
		if subdir not in valid_subdir_options:
			while not cont:
				temp = raw_input('type "no_mask" or "full" or other subdirectory to do the single channel study: ')
				if temp in valid_subdir_options:
					cont = True
					sub_dir = deepcopy(temp)
		CreateDirectoryIfNecessary(self.runlist_dir)
		self.Create_Run_List(sub_dir)
		CreateDirectoryIfNecessary(self.out_dir + '/' + sub_dir + '/' + str(self.run))
		RecreateSoftLink(self.out_dir + '/' + sub_dir + '/' + str(self.run), scratch_path, str(self.run) + '_' + sub_dir)
		self.process_f = subp.Popen(['diamondAnalysis', '-r', '{d}/{sd}/RunList_{r}.ini'.format(d=self.runlist_dir, sd=sub_dir, r=self.run), '-s', self.settings_dir + '/' + sub_dir, '-o', self.out_dir + '/' + sub_dir, '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		while self.process_f.poll() is None:
			pass
		print 'Process finished'
		CloseSubprocess(self.process_f, True, False)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-s', '--settings', dest='settings_f', type='string', help='Settings file containing the information on the run and the analysis to do (e.g. settings.ini)')
	parser.add_option('--first', dest='first', default=False, action='store_true', help='enables first analysis wich has everything un-masked')
	parser.add_option('--normal', dest='normal', default=False, action='store_true', help='enables normal analysis')
	parser.add_option('--singlechannel', dest='singlech', default=False, action='store_true', help='enables single channel study. Requires a preexiting first analysis')

	(options, args) = parser.parse_args()
	settings_f = str(options.settings_f)
	first_ana = bool(options.first)
	normal_ana = bool(options.normal)
	single_ch = bool(options.singlech)

	rd42 = RD42Analysis()
	rd42.ReadInputFile(settings_f)
	rd42.Create_Run_List()
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
