#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
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

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12 - 1}
ped_axis = {'min': 0, 'max': 2**12 - 1}
cm_axis = {'min': -100, 'max': 100}

scratch_path = '/scratch/strip_telescope_tests/runDiego/output'

class RD42Analysis:
	def __init__(self, run, source_dir, numev, runlistdir, settings_dir, pedestal=False, cluster=False, selec=False, autom_all=False):
		print 'Creating RD42Analysis instance for run:', run
		self.run = run
		self.in_dir = source_dir + str(self.run)
		self.out_dir = source_dir + 'output'
		CreateDirectoryIfNecessary(self.out_dir)
		if numev != 0:
			self.num_ev = numev
		else:
			cont = False
			while not cont:
				temp = raw_input('Enter the number of events to analyse: ')
				if IsInt(temp):
					if 0 < temp < 100e6:
						cont = True
				else:
					print 'What you entered is not valid. Try again'
			self.num_ev = deepcopy(temp)
			del temp
		self.runlist_dir = Correct_Path(runlistdir)
		self.settings_dir = Correct_Path(settings_dir)
		self.do_pedestal_ana = pedestal
		self.do_cluster_ana = cluster
		self.do_selec_ana = selec
		self.do_autom_all = autom_all
		self.bar = None
		self.process_f = None
		self.process_e = None
		self.process_o = None
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(999)

	def Create_Run_List_full(self):
		if not os.path.isdir(self.runlist_dir):
			os.makedirs(self.runlist_dir)
		ped = 1 if self.do_pedestal_ana else 0
		clu = 1 if self.do_cluster_ana else 0
		sele = 1 if self.do_selec_ana else 0
		if not os.path.isdir(self.runlist_dir+'/full'):
			os.makedirs(self.runlist_dir+'/full')
		with open(self.runlist_dir + '/full/RunList_'+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t1\t0\t1\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))
		if not os.path.isdir(self.runlist_dir+'/odd'):
			os.makedirs(self.runlist_dir+'/odd')
		with open(self.runlist_dir+'/odd' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))
		if not os.path.isdir(self.runlist_dir+'/even'):
			os.makedirs(self.runlist_dir+'/even')
		with open(self.runlist_dir+'/even' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))

	def Copy_settings_to_even_odd(self):
		self.Check_settings_full()
		if not os.path.isdir(self.settings_dir + '/even'):
			os.makedirs(self.settings_dir + '/even')
		if not os.path.isdir(self.settings_dir + '/odd'):
			os.makedirs(self.settings_dir + '/odd')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/even/')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/odd/')

	def Check_settings_full(self):
		if not os.path.isdir(self.settings_dir + '/full'):
			os.makedirs(self.settings_dir + '/full')
		if not os.path.isfile(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run)):
			CreateDefaultSettingsFile(self.settings_dir + '/full', self.run, self.num_ev)
		Set_Diamond_pattern(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_NC(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_Screened(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_Noisy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Select_fiducial_region(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))

	def Modify_even_odd(self):
		self.Modify_even()
		self.Modify_odd()

	def Modify_even(self):
		Replace_Settings_Line(self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run), 'Dia_channel_screen_channels', 'even')
		# with open(self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run), 'r') as fin:
		# 	with open(self.settings_dir + '/even/settings.new.{r}.ini.'.format(r=self.run), 'w') as fnew:
		# 		for line in fin:
		# 			if not line.startswith('Dia_channel_screen_channels'):
		# 				fnew.write(line)
		# 			else:
		# 				channel_str = line[line.find('{') + 1:line.find('}')].split(',')
		# 				channel_str_new = self.ModifyString(channel_str, 'even')
		# 				fnew.write('Dia_channel_screen_channels = {' + channel_str_new + '}\n')
		# os.remove(self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run))
		# shutil.move(self.settings_dir + '/even/settings.new.{r}.ini.'.format(r=self.run), self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run))

	def Modify_odd(self):
		Replace_Settings_Line(self.settings_dir + '/odd/settings.{r}.ini'.format(r=self.run), 'Dia_channel_screen_channels', 'odd')
		# with open(self.settings_dir + '/odd/settings.{r}.ini'.format(r=self.run), 'r') as fin:
		# 	with open(self.settings_dir + '/odd/settings.new.{r}.ini.'.format(r=self.run), 'w') as fnew:
		# 		for line in fin:
		# 			if not line.startswith('Dia_channel_screen_channels'):
		# 				fnew.write(line)
		# 			else:
		# 				channel_str = line[line.find('{') + 1:line.find('}')].split(',')
		# 				channel_str_new = self.ModifyString(channel_str, 'odd')
		# 				fnew.write('Dia_channel_screen_channels = {' + channel_str_new + '}\n')
		# os.remove(self.settings_dir + '/odd/settings.{r}.ini'.format(r=self.run))
		# shutil.move(self.settings_dir + '/odd/settings.new.{r}.ini.'.format(r=self.run), self.settings_dir + '/odd/settings.{r}.ini'.format(r=self.run))

	# def ModifyString(self, channel_old, type='even'):
	# 	modulo = 1 if type == 'even' else 0
	# 	channel_new = ''
	# 	for i in xrange(1, len(channel_old)):
	# 		if IsInt(channel_old[i - 1]):
	# 			prev = int(channel_old[i - 1])
	# 			if IsInt(channel_old[i]):
	# 				th = int(channel_old[i])
	# 				if th - prev < 2:
	# 					channel_new += str(prev) + ','
	# 				elif prev % 2 == modulo:
	# 					for ch in xrange(prev, th, 2):
	# 						channel_new += str(ch) + ','
	# 				else:
	# 					channel_new += str(prev) + ','
	# 					for ch in xrange(prev + 1, th, 2):
	# 						channel_new += str(ch) + ','
	# 			else:
	# 				temp = channel_old[i].split('-')
	# 				if IsInt(temp[0]):
	# 					th = int(temp[0])
	# 					if th - prev < 2:
	# 						channel_new += str(prev) + ','
	# 					elif prev % 2 == modulo:
	# 						for ch in xrange(prev, th, 2):
	# 							channel_new += str(ch) + ','
	# 					else:
	# 						channel_new += str(prev) + ','
	# 						for ch in xrange(prev + 1, th, 2):
	# 							channel_new += str(ch) + ','
	# 		else:
	# 			channel_new += channel_old[i - 1] + ','
	# 	channel_new = channel_new[:-1]
	# 	return channel_new

	def ModifySettingsOneChannel(self, ch=0):
		CreateDirectoryIfNecessary(self.settings_dir + '/no_mask/channels')
		CreateDirectoryIfNecessary(self.settings_dir + '/no_mask/channels/{c}'.format(c=ch))
		with open(self.settings_dir + '/no_mask/settings.{r}.ini'.format(r=self.run), 'r') as fin:
			with open(self.settings_dir + '/no_mask/channels/{c}/settings.{r}.ini'.format(c=ch, r=self.run), 'w') as fch:
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
		CreateDirectoryIfNecessary(self.runlist_dir)
		if not self.IsMaskRunDone():
			print 'Starting first analysis...'
			with open(self.runlist_dir + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev))
			CreateDefaultSettingsFile(self.settings_dir + '/no_mask', self.run, self.num_ev)
			CreateDirectoryIfNecessary(self.out_dir + '/no_mask')
			RecreateSoftLink(self.out_dir + '/no_mask/' + str(self.run), scratch_path, str(self.run) + '_no_mask')
			self.process_f = subp.Popen(['diamondAnalysis', '-r', '{d}/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/no_mask', '-o', self.out_dir + '/no_mask', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
			while self.process_f.poll() is None:
				continue
			CloseSubprocess(self.process_f, stdin=True, stdout=False)
			print 'Finish first analysis :)'
		if not os.path.isdir(self.out_dir + '/no_mask/{r}/channel_sweep'.format(r=self.run)):
			self.GetIndividualChannelHitmap()
		else:
			skip = True
			for ch in xrange(diaChs):
				if not os.path.isfile(self.out_dir + '/no_mask/{r}/channel_sweep/{c}/{r}/selectionData.{r}.root'.format(r=self.run, c=ch)):
					skip = False
					break
			if not skip:
				self.GetIndividualChannelHitmap()

	def GetIndividualChannelHitmap(self):
		print 'Starting channel occupancy plots...'
		CreateDirectoryIfNecessary(self.out_dir + '/no_mask/{r}/channel_sweep'.format(r=self.run))
		for ch in xrange(diaChs):
			CreateDirectoryIfNecessary(self.out_dir + '/no_mask/{r}/channel_sweep/{c}'.format(c=ch, r=self.run))
			CreateDirectoryIfNecessary(self.out_dir + '/no_mask/{r}/channel_sweep/{c}/{r}'.format(c=ch, r=self.run))
			self.LinkRootFiles(self.out_dir + '/no_mask/' + str(self.run), self.out_dir + '/no_mask/{r}/channel_sweep/{c}/{r}'.format(c=ch, r=self.run), 'cluster')
			self.ModifySettingsOneChannel(ch)
		channel_runs = np.arange(diaChs).reshape([16, 8])
		for bat in xrange(channel_runs.shape[0]):
			procx = []
			for proc in xrange(channel_runs.shape[1]):
				print 'Getting channel:', channel_runs[bat, proc], '...'
				procx.append(subp.Popen(
					['diamondAnalysis', '-r', '{d}/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/no_mask/channels/{c}'.format(c=channel_runs[bat, proc]),
					 '-o', self.out_dir + '/no_mask/{r}/channel_sweep/{c}'.format(c=channel_runs[bat, proc], r=self.run), '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE,
					stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True))
			for proc in xrange(channel_runs.shape[1]):
				while procx[proc].poll() is None:
					# procx[proc].stdout.
					continue
			for proc in xrange(channel_runs.shape[1]):
				CloseSubprocess(procx[proc], stdin=True, stdout=False)
				print 'Done with channel:', channel_runs[bat, proc], ':)'
			del procx

	def Full_Analysis(self):
		CreateDirectoryIfNecessary(self.runlist_dir)
		self.Create_Run_List_full()
		self.Copy_settings_to_even_odd()
		self.Modify_even_odd()
		CreateDirectoryIfNecessary(self.out_dir + '/full/' + str(self.run))
		CreateDirectoryIfNecessary(self.out_dir + '/even/' + str(self.run))
		CreateDirectoryIfNecessary(self.out_dir + '/odd/' + str(self.run))
		RecreateSoftLink(self.out_dir + '/full/' + str(self.run), scratch_path, str(self.run) + '_full')
		RecreateSoftLink(self.out_dir + '/even/' + str(self.run), scratch_path, str(self.run) + '_even')
		RecreateSoftLink(self.out_dir + '/odd/' + str(self.run), scratch_path, str(self.run) + '_odd')
		if os.path.isfile(self.out_dir + '/no_mask/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
			self.LinkRootFiles(self.out_dir + '/no_mask/' + str(self.run), self.out_dir + '/full/' + str(self.run), 'raw')
		self.process_f = subp.Popen(['diamondAnalysis', '-r', '{d}/full/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/full', '-o', self.out_dir + '/full', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		is_finished_cluster = False
		is_even_odd_started = False
		while self.process_f.poll() is None:
			if not is_even_odd_started:
				is_finished_cluster = self.Check_if_clustering_finished()
				if is_finished_cluster:
					self.LinkRootFiles(self.out_dir + '/full/' + str(self.run), self.out_dir + '/even/' + str(self.run), 'cluster', doCopy=False)
					self.LinkRootFiles(self.out_dir + '/full/' + str(self.run), self.out_dir + '/odd/' + str(self.run), 'cluster', doCopy=False)
					self.process_e = subp.Popen(['diamondAnalysis', '-r', '{d}/even/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/even', '-o', self.out_dir + '/even', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True)
					self.process_o = subp.Popen(['diamondAnalysis', '-r', '{d}/odd/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/odd', '-o', self.out_dir + '/odd', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, stdout=open('/dev/null', 'w'), stderr=subp.STDOUT, close_fds=True)
					is_even_odd_started = True
		print 'Full process finished'
		CloseSubprocess(self.process_f, True, False)
		temp_flag = True
		while self.process_e.poll() is None:
			if temp_flag:
				print 'waiting for even events to finish...'
				temp_flag = False
		CloseSubprocess(self.process_e, True, True)
		temp_flag = True
		while self.process_o.poll() is None:
			if temp_flag:
				print 'waiting for even events to finish...'
				temp_flag = False
		CloseSubprocess(self.process_e, True, True)
		print 'Finished :)'

	def Check_if_clustering_finished(self):
		if os.path.isdir(self.out_dir + '/full/{r}/selections'.format(r=self.run)):
			return True
		return False

	def IsMaskRunDone(self):
		no_mask_dir = self.out_dir + '/no_mask/' + str(self.run)
		if os.path.isdir(no_mask_dir):
			if os.path.isfile(no_mask_dir + '/pedestalData.{r}.root'.format(r=self.run)):
				tempf = ro.TFile(no_mask_dir + '/pedestalData.{r}.root'.format(r=self.run), 'read')
				tempt = tempf.Get('pedestalTree')
				if tempt:
					if tempt.GetEntries() >= 10000:
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
			RecreateLink(source + '/rawData.{r}.root'.format(r=self.run), dest, 'rawData.{r}.root'.format(r=self.run))
		if 'pedestal' in cummulative:
			RecreateLink(source + '/pedestalData.{r}.root'.format(r=self.run), dest, 'pedestalData.{r}.root'.format(r=self.run))
			RecreateLink(source + '/pedestalAnalysis', dest, 'pedestalAnalysis')
		if 'cluster' in cummulative:
			RecreateLink(source + '/clusterData.{r}.root'.format(r=self.run), dest, 'clusterData.{r}.root'.format(r=self.run))
			RecreateLink(source + '/etaCorrection.{r}.root'.format(r=self.run), dest, 'etaCorrection.{r}.root'.format(r=self.run))
			RecreateLink(source + '/clustering', dest, 'clustering')
			if os.path.isfile(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run)) and doCopy:
				shutil.copy(source + '/crossTalkCorrectionFactors.{r}.txt'.format(r=self.run), dest + '/')
		if 'selection' in cummulative:
			RecreateLink(source + '/selectionData.{r}.root'.format(r=self.run), dest, 'selectionData.{r}.root'.format(r=self.run))
			RecreateLink(source + '/selectionAnalysis', dest, 'selectionAnalysis')
			RecreateLink(source + '/selections', dest, 'selections')
		if 'align' in cummulative:
			RecreateLink(source + '/alignment.{r}.root'.format(r=self.run), dest, 'alignment.{r}.root'.format(r=self.run))
			RecreateLink(source + '/alignment', dest, 'alignment')
		if 'transparent' in cummulative:
			RecreateLink(source + '/transparentAnalysis', dest, 'transparentAnalysis')
		if os.path.isfile(source + '/index.php'):
			shutil.copyfile(source + '/index.php', dest + '/index.php')
		if os.path.isfile(source + '/overview.html'):
			shutil.copyfile(source + '/overview.html', dest + '/overview.html')
		if os.path.isfile(source + '/results_{r}.res'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.res'.format(r=self.run), dest + '/results_{r}.res'.format(r=self.run))
		if os.path.isfile(source + '/results_{r}.txt'.format(r=self.run)) and doCopy:
			shutil.copyfile(source + '/results_{r}.txt'.format(r=self.run), dest + '/results_{r}.txt'.format(r=self.run))
		RecreateLink(source + '/Results.{r}.root'.format(r=self.run), dest, 'Results.{r}.root'.format(r=self.run))

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
	parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
	parser.add_option('-i', '--input', dest='input', default='.', type='string', help='folder containing the different folder runs')
	# parser.add_option('-o', '--output', dest='output', default='.', type='string', help='foler where to save the structures (will contain, full, even and odd subfolders)')
	parser.add_option('-l', '--runlistdir', dest='runlist', default='~/RunLists', help='folder which contains the RunLists files')
	parser.add_option('-s', '--settingsdir', dest='sett', default='~/settings', help='folder which contains the settings files')
	parser.add_option('-n', '--numevents', dest='numevents', default=0, type='int', help='number of events to analyse')
	parser.add_option('-p', '--pedestal', dest='pedestal', default=False, action='store_true', help='enables pedestal analysis')
	parser.add_option('-c', '--cluster', dest='cluster', default=False, action='store_true', help='enables cluster analysis')
	parser.add_option('-e', '--selection', dest='selection', default=False, action='store_true', help='enables selection analysis')
	parser.add_option('-a', '--all', dest='all', default=False, action='store_true', help='enables all types of analysis. Creates odd, even and full')
	parser.add_option('-x', '--first', dest='first', default=False, action='store_true', help='enables first analysis wich has everything un-masked')

	(options, args) = parser.parse_args()
	run = int(options.run)
	input = str(options.input)
	# output = bool(options.output)
	numev = int(options.numevents)
	runlist = str(options.runlist)
	settings_dir = str(options.sett)
	pedestal = bool(options.pedestal)
	cluster = bool(options.cluster)
	selec = bool(options.selection)
	autom_all = bool(options.all)
	first_ana = bool(options.first)

	rd42 = RD42Analysis(run=run, source_dir=input, numev=numev, runlistdir=runlist, settings_dir=settings_dir, pedestal=pedestal, cluster=cluster, selec=selec, autom_all=autom_all)

	if first_ana:
		rd42.First_Analysis()
	elif autom_all:
		rd42.Full_Analysis()
	else:
		pass # run analysis for full with existing runlist file
	# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

