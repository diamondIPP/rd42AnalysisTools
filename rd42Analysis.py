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

# scratch_path = '/scratch/strip_telescope_tests/runDiego/output'  # at snickers
scratch_path = '/eos/user/d/dsanzbec/scratch/output'  # at lxplus

class RD42Analysis:
	def __init__(self, run, source_dir, output_subdir, numev, runlistdir, settings_dir, pedestal=False, cluster=False, selec=False, doAlign=False, transparent=False, evini=0, numEvsAna=0, deleteold=False):
		print 'Creating RD42Analysis instance for run:', run
		self.run = run
		self.in_dir = source_dir + str(self.run)
		self.out_dir = source_dir + 'output'
		self.subdir = output_subdir
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
		self.firstev = evini
		self.num_ev_ana = self.num_ev - self.firstev if numEvsAna == 0 else numEvsAna
		self.runlist_dir = Correct_Path(runlistdir)
		self.settings_dir = Correct_Path(settings_dir)
		self.do_pedestal_ana = pedestal
		self.do_cluster_ana = cluster
		self.do_selec_ana = selec
		self.do_align = doAlign
		self.do_trans = transparent
		self.bar = None
		self.process_f = None
		self.process_e = None
		self.process_o = None
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(999)
		self.deleteold = deleteold

	def Create_Run_List(self, subdir='full'):
		if not os.path.isdir(self.runlist_dir):
			os.makedirs(self.runlist_dir)
		ped = 1 if self.do_pedestal_ana else 0
		clu = 1 if self.do_cluster_ana else 0
		sele = 1 if self.do_selec_ana else 0
		alig = 1 if self.do_align else 0
		tran = 1 if self.do_trans else 0
		if not os.path.isdir(self.runlist_dir+'/{f}'.format(f=subdir)):
			os.makedirs(self.runlist_dir+'/{f}'.format(f=subdir))
		with open(self.runlist_dir + '/{f}/RunList_'.format(f=subdir)+str(self.run)+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t{al}\t0\t{t}\n#\n'.format(r=self.run, n=self.num_ev_ana, p=ped, c=clu, s=sele, al=alig, t=tran))
		if subdir == 'full':
			if not os.path.isdir(self.runlist_dir+'/odd'):
				os.makedirs(self.runlist_dir+'/odd')
			with open(self.runlist_dir+'/odd' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.num_ev_ana, p=ped, c=clu, s=sele))
			if not os.path.isdir(self.runlist_dir+'/even'):
				os.makedirs(self.runlist_dir+'/even')
			with open(self.runlist_dir+'/even' + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n#\n'.format(r=self.run, n=self.num_ev_ana, p=ped, c=clu, s=sele))

	def Copy_settings_to_even_odd(self):
		self.Check_settings_full()
		if not os.path.isdir(self.settings_dir + '/even'):
			os.makedirs(self.settings_dir + '/even')
		if not os.path.isdir(self.settings_dir + '/odd'):
			os.makedirs(self.settings_dir + '/odd')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/even/')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/odd/')

	def Check_settings_full(self):
		print 'Checking settings...\n'
		if not os.path.isdir(self.settings_dir + '/full'):
			os.makedirs(self.settings_dir + '/full')
		if not os.path.isfile(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run)):
			CreateDefaultSettingsFile(self.settings_dir + '/full', self.run, self.num_ev, ev_ini=self.firstev, num_evs_ana=self.num_ev_ana)
		Set_voltage(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Set_Diamond_pattern(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Set_Diamond_name(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_NC(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_Screened(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Mark_channels_as_Noisy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		if self.do_align:
			Mark_channels_as_no_alignment(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))
		Select_fiducial_region(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run))

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
		CreateDirectoryIfNecessary(self.runlist_dir)
		if not self.IsMaskRunDone():
			print 'Starting first analysis...'
			with open(self.runlist_dir + '/RunList_'+str(self.run)+'.ini', 'w') as rlf:
				rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev))
			if not os.path.isfile(self.settings_dir + '/no_mask/settings.' + str(self.run) + '.ini'):
				CreateDefaultSettingsFile(self.settings_dir + '/no_mask', self.run, self.num_ev)
			CreateDirectoryIfNecessary(self.out_dir + '/no_mask')
			CreateDirectoryIfNecessary(self.out_dir + '/no_mask/' + str(self.run))
			RecreateSoftLink(self.out_dir + '/no_mask/' + str(self.run), scratch_path, str(self.run) + '_no_mask')
			self.process_f = subp.Popen(['diamondAnalysis', '-r', '{d}/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/no_mask', '-o', self.out_dir + '/no_mask', '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
			while self.process_f.poll() is None:
				continue
			CloseSubprocess(self.process_f, stdin=True, stdout=False)
			print 'Finish first analysis :)'
		else:
			print 'First analysis already exists :)'


	def GetIndividualChannelHitmap(self, subdir=''):
		print 'Starting channel occupancy plots...'
		cont = False
		sub_dir = subdir
		valid_subdir_options = ['no_mask', 'full']
		if subdir not in valid_subdir_options:
			while not cont:
				temp = raw_input('type "no_mask" or "full" or other subdirectory to do the single channel study: ')
				if temp in valid_subdir_options:
					cont = True
					sub_dir = deepcopy(temp)

		CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep'.format(sd=sub_dir, r=self.run))
		for ch in xrange(diaChs):
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=sub_dir, c=ch, r=self.run))
			CreateDirectoryIfNecessary(self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=sub_dir, c=ch, r=self.run))
			self.LinkRootFiles(self.out_dir + '/{sd}/' + str(self.run), self.out_dir + '/{sd}/{r}/channel_sweep/{c}/{r}'.format(sd=sub_dir, c=ch, r=self.run), 'cluster')
			self.ModifySettingsOneChannel(ch, sub_dir)
		num_cores = mp.cpu_count()
		num_parrallel = 1 if num_cores < 4 else 2 if num_cores < 8 else 4 if num_cores < 16 else 8 if num_cores < 32 else 16
		num_jobs = diaChs / num_parrallel
		channel_runs = np.arange(diaChs).reshape([num_jobs, num_parrallel])
		for bat in xrange(channel_runs.shape[0]):
			procx = []
			for proc in xrange(channel_runs.shape[1]):
				print 'Getting channel:', channel_runs[bat, proc], '...'
				procx.append(subp.Popen(
					['diamondAnalysis', '-r', '{d}/RunList_{r}.ini'.format(sd=sub_dir, d=self.runlist_dir, r=self.run), '-s', self.settings_dir + '/{sd}/channels/{c}'.format(sd=sub_dir, c=channel_runs[bat, proc]),
					 '-o', self.out_dir + '/{sd}/{r}/channel_sweep/{c}'.format(sd=sub_dir, c=channel_runs[bat, proc], r=self.run), '-i', self.in_dir], bufsize=-1, stdin=subp.PIPE,
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
		self.Create_Run_List('full')
		self.Copy_settings_to_even_odd()
		self.Modify_even_odd()
		if self.deleteold:
			DeleteDirectoryContents(self.out_dir + '/full/' + str(self.run))
			DeleteDirectoryContents(self.out_dir + '/even/' + str(self.run))
			DeleteDirectoryContents(self.out_dir + '/odd/' + str(self.run))
		CreateDirectoryIfNecessary(self.out_dir + '/full/' + str(self.run))
		CreateDirectoryIfNecessary(self.out_dir + '/even/' + str(self.run))
		CreateDirectoryIfNecessary(self.out_dir + '/odd/' + str(self.run))
		RecreateSoftLink(self.out_dir + '/full/' + str(self.run), scratch_path, str(self.run) + '_full')
		RecreateSoftLink(self.out_dir + '/even/' + str(self.run), scratch_path, str(self.run) + '_even')
		RecreateSoftLink(self.out_dir + '/odd/' + str(self.run), scratch_path, str(self.run) + '_odd')
		if os.path.isfile(self.out_dir + '/no_mask/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
			if self.firstev == 0 and self.num_ev_ana == self.num_ev:
				print 'Linking raw file from no_mask...',; sys.stdout.flush()
				self.LinkRootFiles(self.out_dir + '/no_mask/' + str(self.run), self.out_dir + '/full/' + str(self.run), 'raw')
				print 'Done'
			else:
				recreate = True
				if os.path.isfile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)) or os.path.islink(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run)):
					tempf0 = ro.TFile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'READ')
					tempt0 = tempf0.Get('rawTree')
					recreate = False if tempt0.GetEntries() == self.num_ev_ana else True
					tempf0.Close()
				if recreate:
					tempf = ro.TFile(self.out_dir + '/no_mask/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'READ')
					tempt = tempf.Get('rawTree')
					print 'Extracting only {eva} events starting from {evi} for analysis...'.format(eva=self.num_ev_ana, evi=self.firstev),; sys.stdout.flush()
					leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.firstev, eva=self.num_ev_ana))
					# leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{evf}-{evi})<=({evf}-{evi})'.format(evi=self.firstev, evf=self.firstev+self.num_ev_ana-1))
					while leng > tempt.GetEstimate():
						tempt.SetEstimate(leng)
						leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.firstev, eva=self.num_ev_ana))
					# leng = tempt.Draw('>>evlist', 'abs(2*EventNumber-{evf}-{evi})<=({evf}-{evi})'.format(evi=self.firstev, evf=self.firstev+self.num_ev_ana-1))
					evlist = ro.gDirectory.Get('evlist')
					tempt.SetEventList(evlist)
					tempnf = ro.TFile(self.out_dir + '/full/' + str(self.run) + '/rawData.{r}.root'.format(r=self.run), 'RECREATE')
					tempnt = tempt.CopyTree('')
					tempnt.Write()
					tempnf.Close()
					tempf.Close()
					print 'Done'
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
	parser.add_option('-r', '--run', dest='run', default=0, type='int', help='Run to be analysed (e.g. 22011)')
	parser.add_option('-i', '--input', dest='input', default='.', type='string', help='folder containing the different folder runs')
	parser.add_option('-o', '--output', dest='output', default='', type='string', help='subfoler where to save the structures (e.g. full, no_mask, poly, etc...)')
	parser.add_option('-l', '--runlistdir', dest='runlist', default='~/RunLists', help='folder which contains the RunLists files')
	parser.add_option('-s', '--settingsdir', dest='sett', default='~/settings', help='folder which contains the settings files')
	parser.add_option('-n', '--numevents', dest='numevents', default=0, type='int', help='number of events to analyse')
	parser.add_option('--pedestal', dest='pedestal', default=False, action='store_true', help='enables pedestal analysis')
	parser.add_option('--cluster', dest='cluster', default=False, action='store_true', help='enables cluster analysis')
	parser.add_option('--selection', dest='selection', default=False, action='store_true', help='enables selection analysis')
	parser.add_option('--alignment', dest='alignment', default=False, action='store_true', help='enables alignment')
	parser.add_option('--transparent', dest='transparent', default=False, action='store_true', help='enables transparent analysis')
	parser.add_option('--full', dest='full', default=False, action='store_true', help='enables all types of analysis. Creates odd, even and full')
	parser.add_option('--first', dest='first', default=False, action='store_true', help='enables first analysis wich has everything un-masked')
	parser.add_option('-x', '--singlechannel', dest='singlech', default=False, action='store_true', help='enables single channel study. Requires a preexiting first analysis')
	parser.add_option('--firstevent', dest='firstevent', default=0, type='int', help='first event to analyse')
	parser.add_option('--numevsana', dest='numevsana', default=0, type='int', help='number of events to analyse')
	parser.add_option('--deleteold', dest='deleteold', default=False, action='store_true', help='deletes previous analysis and creates a new one')

	(options, args) = parser.parse_args()
	run = int(options.run)
	input = str(options.input)
	output = bool(options.output)
	numev = int(options.numevents)
	runlist = str(options.runlist)
	settings_dir = str(options.sett)
	pedestal = bool(options.pedestal)
	cluster = bool(options.cluster)
	selec = bool(options.selection)
	alig = bool(options.alignment)
	tran = bool(options.transparent)
	fullana = bool(options.full)
	first_ana = bool(options.first)
	single_ch = bool(options.singlech)
	firstev = int(options.firstevent)
	numevsana = int(options.numevsana)
	deleteold = bool(options.deleteold)

	rd42 = RD42Analysis(run=run, source_dir=input, output_subdir=output, numev=numev, runlistdir=runlist, settings_dir=settings_dir, pedestal=pedestal, cluster=cluster, selec=selec, doAlign=alig, transparent=tran, evini=firstev, numEvsAna=numevsana, deleteold=deleteold)

	if first_ana:
		print 'Starting first analysis (no_mask)...\n'
		rd42.First_Analysis()
	elif fullana:
		print 'Starting full analysis (full)...\n'
		rd42.Full_Analysis()
	else:
		pass  # run analysis for full with existing runlist file
	if single_ch:
		if first_ana:
			rd42.GetIndividualChannelHitmap(subdir='no_mask')
		elif fullana:
			rd42.GetIndividualChannelHitmap(subdir='full')
		else:
			rd42.GetIndividualChannelHitmap()
# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

