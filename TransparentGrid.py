#!/usr/bin/env python
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
from copy import deepcopy
from collections import OrderedDict
import os, sys, shutil
from Utils import *
import cPickle as pickle
from Langaus import LanGaus
from GridAreas import GridAreas

ph_bins_options = np.array((1, 2, 4, 5, 10, 16, 20, 25, 32, 40, 50, 80, 100, 125, 160, 200, 250, 400, 500, 800, 1000, 2000), 'uint16')

class TransparentGrid:
	def __init__(self, dir='', run=25209, cellsize=50.0):
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(99)
		ro.gStyle.SetOptStat(11)
		ro.gStyle.SetOptFit(111)
		ro.gStyle.SetPaintTextFormat(".0f")
		# ro.TFormula.SetMaxima(100000)
		self.run = run
		self.dir = os.path.abspath(os.path.expanduser(os.path.expandvars(dir)))
		self.trans_file = None
		self.trans_tree = None
		self.align_file = None
		self.align_obj = None
		self.pkl = None
		self.pkl_sbdir = 'default'
		self.loaded_pickle = False
		self.align_info = {'xoff': float(0), 'phi': float(0)}
		self.cluster_size, self.num_strips = 0, 0
		self.num_cols = 19 if cellsize == 50 else 13
		self.ch_ini = 0
		self.ch_end = 84
		self.phbins = 200
		self.phmin = 0
		self.phmax = 4000
		self.phbins_neg = 200
		self.phmin_neg = -2000
		self.phmax_neg = 2000
		self.neg_cut = 5
		self.col_pitch = cellsize
		self.cell_resolution = 50.0 / 25 if self.col_pitch == 50 else 100.0 / 51
		self.delta_offset_threshold = 0.01  # in mum
		# self.saturated_ADC = 4095
		if self.run in [25207, 25208, 25209, 25210]:
			if self.run in [25207, 25208, 25210]: print 'Using settings for run 25209. Results might not be correctly aligned'
			self.row_info_diamond = {'num': 27, 'pitch': 50.0, 'x_off': 0.5065, 'y_off': 18.95, '0': 3117.70005, 'up': 4467.70005} if self.col_pitch == 50 else {'num': 24, 'pitch': 100.0, 'x_off': 0.4976, 'y_off': 52.6, '0': 3049.6, 'up': 5449.6}
		elif self.run in [25202, 25203, 25204, 25205, 25206]:
			if self.run in [25202, 25203, 25206]: print 'Using settings for run 25205. Results might not be correctly aligned'
			self.row_info_diamond = {'num': 27, 'pitch': 50.0, 'x_off': 0.4922, 'y_off': 47.978033, '0': 3148.785, 'up': 4498.785} if self.col_pitch == 50 else {'num': 24, 'pitch': 100.0, 'x_off': 0.4976, 'y_off': 79.803, '0': 3076.115, 'up': 5476.115}
		self.bins_per_ch_x = 3 if self.col_pitch == 50 else 5
		self.bins_per_ch_y = 3 if self.col_pitch == 50 else 5
		self.length_central_region = 30 if self.col_pitch == 50 else 40
		self.conv_steps = 1000
		self.sigma_conv = 5
		self.vertical_lines_diamond = []
		self.vertical_lines_diamond_tline = []
		self.horizontal_lines_diamond = []
		self.horizontal_lines_diamond_tline = []
		self.mpshift = -0.22278298
		self.canvas = {}
		self.line = {}
		self.profile = {}
		self.histo = {}
		self.names = []
		self.tcutgs_diamond = {}
		self.tcutg_diamond_center = None
		self.tcutgs_diamond_center = {}
		self.gridAreas = None
		self.temph = None
		self.langaus = {}
		self.gridTextDiamond = None
		self.doubleLangaus = {}
		self.lg1, self.lg2 = ro.TF1(), ro.TF1()
		self.CheckFoldersAndFiles()
		self.OpenFileAndGetTree()
		self.FindDiamondChannelLimits()

	def CheckFoldersAndFiles(self):
		if not os.path.isdir(self.dir):
			ExitMessage('The path to the directory "{d}" does not exist. Exiting...'.format(d=self.dir), code=os.EX_NOINPUT)
		if not os.path.isdir('{d}/{r}'.format(d=self.dir, r=self.run)):
			ExitMessage('There is no run {r} in the directory "{d}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)
		if not os.path.isfile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run)):
			ExitMessage('There is no transparent root file "transparent.{r}.root" in the directory "{d}/{r}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)

	def OpenFileAndGetTree(self):
		self.trans_file = ro.TFile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run), 'r')
		self.trans_tree = self.trans_file.Get('transparentTree')
		tempnh = ro.TH1F('tempn', 'tempn', 128, 0, 128)
		tempch = ro.TH1F('tempc', 'tempc', 128, 0, 128)
		if self.trans_tree.FindLeaf('numStrips'):
			self.trans_tree.Draw('numStrips>>tempn','transparentEvent', 'goff')
			self.num_strips = int(np.floor(tempnh.GetMean() + 0.5))
		if self.trans_tree.FindLeaf('clusterSize'):
			self.trans_tree.Draw('numStrips>>tempc', 'transparentEvent', 'goff')
			self.cluster_size = int(np.floor(tempch.GetMean() + 0.5))

	def FindDiamondChannelLimits(self):
		temph = ro.TH1F('temph', 'temph', 128, -0.5,127.5)
		self.trans_tree.Draw('diaChXPred>>temph', 'transparentEvent', 'goff')
		self.ch_ini = int(temph.GetBinCenter(temph.FindFirstBinAbove(1, 1)))
		self.ch_end = int(temph.GetBinCenter(temph.FindLastBinAbove(1, 1)))
		self.num_cols = self.ch_end - self.ch_ini + 1

	def SavePickle(self):
		object_dic = {}
		object_dic['row_info_diamond'] = self.row_info_diamond
		object_dic['align_info'] = self.align_info
		object_dic['num_cols'] = self.num_cols
		object_dic['ch_ini'] = self.ch_ini
		object_dic['ch_end'] = self.ch_end
		object_dic['phbins'] = self.phbins
		object_dic['phmin'] = self.phmin
		object_dic['phmax'] = self.phmax
		object_dic['phbins_neg'] = self.phbins_neg
		object_dic['phmin_neg'] = self.phmin_neg
		object_dic['phmax_neg'] = self.phmax_neg
		object_dic['neg_cut'] = self.neg_cut
		object_dic['col_pitch'] = self.col_pitch
		object_dic['cell_resolution'] = self.cell_resolution
		# object_dic['saturated_ADC'] = self.saturated_ADC
		object_dic['bins_per_ch_x'] = self.bins_per_ch_x
		object_dic['bins_per_ch_y'] = self.bins_per_ch_y
		object_dic['length_central_region'] = self.length_central_region
		object_dic['conv_steps'] = self.conv_steps
		object_dic['sigma_conv'] = self.sigma_conv
		object_dic['mpshift'] = self.mpshift

		if not os.path.isdir('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
		picklepath = '{d}/{r}/{s}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)

		pickle.dump(object_dic, open(picklepath, 'wb'))

		if not os.path.isdir('{d}/{r}/default'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/default'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
			picklepath = '{d}/{r}/default/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
			pickle.dump(object_dic, open(picklepath, 'wb'))

		print 'Saved pickle :D'

	def LoadPickle(self):
		picklepath = '{d}/{r}/{s}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
		picklepath2 = '{d}/{r}/default/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
		if not os.path.isfile(picklepath):
			picklepath = picklepath2
			print 'Trying to load pickle inside default...'
		if os.path.isfile(picklepath):
			with open(picklepath, 'rb') as pkl:
				self.pkl = pickle.load(pkl)
				self.loaded_pickle = (picklepath != picklepath2)
				if self.loaded_pickle:
					print 'Loaded from pickle in {sd} :D'.format(sd=self.pkl_sbdir)
				else:
					print 'Loaded from pickle in default :D'

	def UnfoldPickle(self):
		if 'row_info_diamond' in self.pkl.keys():
			self.row_info_diamond = self.pkl['row_info_diamond']
		if 'align_info' in self.pkl.keys():
			self.align_info = self.pkl['align_info']
		if 'num_cols' in self.pkl.keys():
			self.num_cols = self.pkl['num_cols']
		if 'ch_ini' in self.pkl.keys():
			self.ch_ini = self.pkl['ch_ini']
		if 'ch_end' in self.pkl.keys():
			self.ch_end = self.pkl['ch_end']
		if 'phbins' in self.pkl.keys():
			self.phbins = self.pkl['phbins']
		if 'phmin' in self.pkl.keys():
			self.phmin = self.pkl['phmin']
		if 'phmax' in self.pkl.keys():
			self.phmax = self.pkl['phmax']
		if 'phbins_neg' in self.pkl.keys():
			self.phbins_neg = self.pkl['phbins_neg']
		if 'phmin_neg' in self.pkl.keys():
			self.phmin_neg = self.pkl['phmin_neg']
		if 'phmax_neg' in self.pkl.keys():
			self.phmax_neg = self.pkl['phmax_neg']
		if 'neg_cut' in self.pkl.keys():
			self.neg_cut = self.pkl['neg_cut']
		if 'col_pitch' in self.pkl.keys():
			self.col_pitch = self.pkl['col_pitch']
		if 'cell_resolution' in self.pkl.keys():
			self.cell_resolution = self.pkl['cell_resolution']
		# if 'saturated_ADC' in self.pkl.keys():
		# 	self.saturated_ADC = self.pkl['saturated_ADC']
		if 'bins_per_ch_x' in self.pkl.keys():
			self.bins_per_ch_x = self.pkl['bins_per_ch_x']
		if 'bins_per_ch_y' in self.pkl.keys():
			self.bins_per_ch_y = self.pkl['bins_per_ch_y']
		if 'length_central_region' in self.pkl.keys():
			self.length_central_region = self.pkl['length_central_region']
		if 'conv_steps' in self.pkl.keys():
			self.conv_steps = self.pkl['conv_steps']
		if 'sigma_conv' in self.pkl.keys():
			self.sigma_conv = self.pkl['sigma_conv']
		if 'mpshift' in self.pkl.keys():
			self.mpshift = self.pkl['mpshift']

	def SetLines(self, try_align=True):
		self.LoadPickle()
		if self.loaded_pickle:
			self.UnfoldPickle()
		elif try_align:
			self.FindHorizontalParametersThroughAlignment()
		else:
			self.AskUserLowerYLines()
			self.CreateLines()
		self.gridAreas = GridAreas(self.num_cols, self.row_info_diamond['num'])

	def FindHorizontalParametersThroughAlignment(self):
		self.LoadAlignmentParameters()

	def LoadAlignmentParameters(self):
		if os.path.isfile('{d}/{r}/alignment.{r}.root'.format(d=self.dir, r=self.run)):
			self.align_file = ro.TFile('{d}/{r}/alignment.{r}.root'.format(d=self.dir, r=self.run), 'r')
			self.align_obj = self.align_file.Get('alignment')
			self.align_info['xoff'] = self.align_obj.GetXOffset(4)
			self.align_info['phi'] = self.align_obj.GetPhiXOffset(4)

	def FindXandYOffests(self):
		if not self.loaded_pickle:
			self.FindUpperAndLowerLines()
			self.FindBinningAndResolution()
			self.row_info_diamond['y_off'] += 3.0 * self.row_info_diamond['pitch'] / 2.0
			delta_y = self.row_info_diamond['pitch']
			proj_width = self.col_pitch / 5.0  # in mum
			proj_bins = int(np.floor(proj_width / self.cell_resolution + 0.5, dtype='float64'))
			proj_bins = proj_bins if proj_bins % 2 == 1 else proj_bins + 1
			proj_low = int(np.floor(np.floor(self.col_pitch / self.cell_resolution + 0.5, dtype='float64') / 2.0 + 0.5, dtype='float64') - (np.floor(proj_bins / 2.0 + 0.5, dtype='float64') - 1) + 1)
			proj_high = proj_low + proj_bins - 1
			iteration = 0
			min_delta = self.row_info_diamond['pitch']
			y_off_shifted, y_min = 0.0, 0.0
			while abs(delta_y) > self.delta_offset_threshold and iteration < 50:
				self.row_info_diamond['y_off'] -= delta_y * (np.exp(-iteration) * 0.9 + 0.1)
				self.DrawProfile2DDiamondCellOverlay('y_off_alignment', 'clusterCharge1', 'good')
				h_proj_y = self.profile['y_off_alignment'].ProjectionY('y_off_alignment_py', proj_low, proj_high, 'e')
				h_proj_y.GetXaxis().SetRangeUser(0, self.row_info_diamond['pitch'])
				miny = h_proj_y.GetBinCenter(h_proj_y.GetMinimumBin())
				miny = miny if abs(miny - self.col_pitch / 2.0) > self.cell_resolution * 2.5 else self.col_pitch / 2.0
				fit_py = h_proj_y.Fit('pol2', 'QEFSN', '', max(miny - 2.5 * self.cell_resolution - 1e-12, 0), min(miny + 2.5 * self.cell_resolution + 1e-12, self.row_info_diamond['pitch']))
				y_min = -fit_py.Parameter(1) / (2 * fit_py.Parameter(2))
				delta_y = self.row_info_diamond['pitch'] / 2.0 - y_min
				if abs(delta_y) < min_delta:
					min_delta = abs(delta_y)
					y_off_shifted = self.row_info_diamond['y_off']
				iteration += 1
				print iteration, y_min, delta_y
			print 'final', min_delta
			self.row_info_diamond['y_off'] = y_off_shifted - self.row_info_diamond['pitch'] / 2.0

			self.row_info_diamond['x_off'] += 3.0 / 2.0
			delta_x = self.col_pitch
			proj_width = self.row_info_diamond['pitch'] / 5.0  # in mum
			proj_bins = int(np.floor(proj_width / self.cell_resolution + 0.5, dtype='float64'))
			proj_bins = proj_bins if proj_bins % 2 == 1 else proj_bins + 1
			proj_low = int(np.floor(np.floor(self.row_info_diamond['pitch'] / self.cell_resolution + 0.5, dtype='float64') / 2.0 + 0.5, dtype='float64') - (np.floor(proj_bins / 2.0 + 0.5, dtype='float64') - 1) + 1)
			proj_high = proj_low + proj_bins - 1
			iteration = 0
			min_delta = self.col_pitch
			x_off_shifted, x_min = 0.0, 0.0
			while abs(delta_x) > self.delta_offset_threshold and iteration < 50:
				self.row_info_diamond['x_off'] -= np.divide(delta_x, self.col_pitch, dtype='float64') * (np.exp(-iteration) * 0.9 + 0.1)
				self.DrawProfile2DDiamondCellOverlay('x_off_alignment', 'clusterCharge1', 'good')
				h_proj_x = self.profile['x_off_alignment'].ProjectionX('x_off_alignment_px', proj_low, proj_high)
				h_proj_x.GetXaxis().SetRangeUser(0, self.col_pitch)
				minx = h_proj_x.GetBinCenter(h_proj_x.GetMinimumBin())
				minx = minx if abs(minx - self.row_info_diamond['pitch'] / 2.0) > self.cell_resolution * 2.5 else self.row_info_diamond['pitch'] / 2.0
				fit_px = h_proj_x.Fit('pol2', 'QEFSN', '', max(minx - 2.5 * self.cell_resolution - 1e-12, 0), min(minx + 2.5 * self.cell_resolution + 1e-12, self.col_pitch))
				x_min = -fit_px.Parameter(1) / (2 * fit_px.Parameter(2))
				delta_x = self.col_pitch / 2.0 - x_min
				if abs(delta_x) < min_delta:
					min_delta = abs(delta_x)
					x_off_shifted = self.row_info_diamond['x_off']
				iteration += 1
				print iteration, x_min, delta_x
			print 'final', min_delta
			self.row_info_diamond['x_off'] = x_off_shifted - 0.5
			ipdb.set_trace()

			self.SavePickle()

	def FindUpperAndLowerLines(self):
		self.DrawProfile2DDiamond('vertical_limits_profile', 'clusterChargeN', '', True)
		xbinmin, xbinmax = self.profile['vertical_limits_profile'].GetXaxis().FindBin(self.ch_ini - 0.5), self.profile['vertical_limits_profile'].GetXaxis().FindBin(self.ch_ini - 0.5) + self.num_cols * self.bins_per_ch_x - 1
		prof_proj_y = self.profile['vertical_limits_profile'].ProjectionY('vertical_limits_profile_py', xbinmin, xbinmax, 'e')
		minbiny, maxbiny = prof_proj_y.FindFirstBinAbove(), prof_proj_y.FindLastBinAbove()
		for biny in xrange(maxbiny, int(prof_proj_y.GetXaxis().GetNbins())):
			if prof_proj_y.GetBinContent(biny) != 0:
				maxbiny = biny
		miny, maxy = prof_proj_y.GetXaxis().GetBinLowEdge(minbiny), prof_proj_y.GetXaxis().GetBinLowEdge(maxbiny + 1)
		func = ro.TF1('box_fcn', '[0]*(TMath::Erf((x-[1])/[2])+1)/2-[3]*(TMath::Erf((x-[4])/[5])+1)/2+[6]', miny, maxy)
		func.SetNpx(int(self.row_info_diamond['num'] * self.bins_per_ch_y * 10))
		zmin, zmax = prof_proj_y.GetMinimum(), prof_proj_y.GetMaximum()
		y1bin, y2bin = prof_proj_y.FindFirstBinAbove((zmin + zmax) / 2.0), prof_proj_y.FindLastBinAbove((zmin + zmax) / 2.0) + 1
		y1, y2 = prof_proj_y.GetXaxis().GetBinCenter(y1bin), prof_proj_y.GetXaxis().GetBinCenter(y2bin)
		z0, z1, z2 = prof_proj_y.GetBinContent(int((minbiny + y1bin) / 2.0)), prof_proj_y.GetBinContent(int((y1bin + y2bin) / 2.0)), prof_proj_y.GetBinContent(int((maxbiny + y2bin) / 2.0))
		func.SetParLimits(0, abs(z1 - z0) / 10.0, 2.0 * abs(z1 - z0))
		func.SetParLimits(1, y1 - 200, y1 + 200)
		func.SetParLimits(2, 0.1, 200)
		func.SetParLimits(3, abs(z1 - z2) / 10.0, 2.0 * abs(z1 - z2))
		func.SetParLimits(4, y2 - 200, y2 + 200)
		func.SetParLimits(5, 0.1, 200)
		func.SetParLimits(6, -2.0 * abs(z0), 10 * abs(z0))
		params = np.array((abs(z1 - z0), y1, 20, abs(z1 - z2), y2, 20, z0), 'float64')
		func.SetParameters(params)
		fit_prof_proj_y = prof_proj_y.Fit('box_fcn', 'QEBMS', 'goff', prof_proj_y.GetBinLowEdge(int((minbiny + y1bin) / 2.0)), prof_proj_y.GetBinLowEdge(int((maxbiny + y2bin) / 2.0)))
		params = np.array((fit_prof_proj_y.Parameter(0), fit_prof_proj_y.Parameter(1), fit_prof_proj_y.Parameter(2), fit_prof_proj_y.Parameter(3), fit_prof_proj_y.Parameter(4), fit_prof_proj_y.Parameter(5), fit_prof_proj_y.Parameter(6)), 'float64')
		func.SetParameters(params)
		fit_prof_proj_y = prof_proj_y.Fit('box_fcn', 'QEBMS', 'goff', (miny + fit_prof_proj_y.Parameter(1))/2.0, (maxy + fit_prof_proj_y.Parameter(4))/2.0)
		extra_y_ch_sharing = fit_prof_proj_y.Parameter(4) - fit_prof_proj_y.Parameter(1) - self.row_info_diamond['num'] * self.row_info_diamond['pitch']
		if extra_y_ch_sharing < 0:
			ExitMessage('Obtained negative vertical extents due to charge sharing ({v}). Check number of rows, or pitch. Exiting'.format(v=extra_y_ch_sharing), os.EX_DATAERR)
		self.row_info_diamond['0'] = fit_prof_proj_y.Parameter(1) + extra_y_ch_sharing / 2.0
		self.row_info_diamond['up'] = fit_prof_proj_y.Parameter(4) - extra_y_ch_sharing / 2.0

	def FindBinningAndResolution(self):
		if self.gridAreas:
			if len(self.gridAreas.goodAreas_index) == 0:
				self.SelectGoodAndBadByThreshold()
		else:
			self.SetLines()
			self.CreateTCutGs()
			self.SelectGoodAndBadByThreshold()
		self.DrawPHGoodAreas('binning_temp', 'clusterCharge1')
		histo_entries = float(self.histo['binning_temp'].GetEntries())
		temp = np.abs(np.subtract(ph_bins_options, histo_entries, dtype='float32'), dtype='float32')
		self.phbins = ph_bins_options[temp.argmin()]
		self.phbins_neg = ph_bins_options[temp.argmin()]
		cell_bins = int(np.floor(np.sqrt(histo_entries * ((self.col_pitch / 2.0) ** 2) / 4500.0, dtype='float64') + 0.5))
		self.cell_resolution = np.divide(self.col_pitch, cell_bins, dtype='float64') if cell_bins % 2 == 1 else np.divide(50.0, cell_bins + 1, dtype='float64')

	def AskUserLowerYLines(self):
		do_diamond = raw_input('Enter 1 if you want to enter the lower y line parameters of the plots in diamond space')
		if bool(int(do_diamond)):
			self.AskUserDiamondLineParameters()

	def AskUserDiamondLineParameters(self):
		self.row_info_diamond['0'] = self.GetFromUser('Enter the y axis intercept in silicon space for the lower detector limit (scalar between 0 and 12800 in um): ', typ='float', limits=[0, 12800])
		self.row_info_diamond['pitch'] = self.GetFromUser('Enter the effective vertical pitch in um: ', typ='float', limits=[0, 20000])
		self.row_info_diamond['x_off'] = self.GetFromUser('Enter the offset for dia X ch for overlay plots (scalar between -1 and 1): ', typ='float', limits=[-1, 1])
		self.row_info_diamond['y_off'] = self.GetFromUser('Enter the offset for dia Y in um for overlay plots (scalar between -{p} and {p}): '.format(p=self.row_info_diamond['pitch']), typ='float', limits=[-self.row_info_diamond['pitch'], self.row_info_diamond['pitch']])
		self.row_info_diamond['num'] = self.GetFromUser('Enter the number of rows: ', typ='int', limits=[1, 1000])

	def GetFromUser(self, message, typ, limits=[]):
		cont = False
		tempv = 0
		while not cont:
			tempv = raw_input(message)
			if typ == 'int':
				if IsInt(tempv):
					tempv = int(tempv)
					if len(limits) == 2:
						if limits[0] <= tempv <= limits[1]:
							cont = True
					else:
						cont = True
			if typ == 'float':
				if IsFloat(tempv):
					tempv = float(tempv)
					if len(limits) == 2:
						if limits[0] <= tempv <= limits[1]:
							cont = True
					else:
						cont = True
		return tempv

	def CreateLines(self):
		linev = self.GetVerticalLineDiamond(x=self.ch_ini - 0.5)
		lineh = self.GetHorizontalLineDiamond(y=self.row_info_diamond['0'])
		self.vertical_lines_diamond.append(linev)
		self.horizontal_lines_diamond.append(lineh)
		for col in xrange(self.num_cols):
			linev = self.GetVerticalLineDiamond(x=self.ch_ini + col + 0.5)
			self.vertical_lines_diamond.append(linev)
		for row in xrange(self.row_info_diamond['num']):
			lineh = self.GetHorizontalLineDiamond(y=self.row_info_diamond['0'] + (row + 1) * self.row_info_diamond['pitch'])
			self.horizontal_lines_diamond.append(lineh)

	def GetVerticalLineDiamond(self, x):
		return {0: {'x': x, 'y': self.row_info_diamond['0']}, 1: {'x': x, 'y': self.row_info_diamond['0'] + self.row_info_diamond['num'] * self.row_info_diamond['pitch']}}

	def GetHorizontalLineDiamond(self, y):
		return {0: {'x': self.ch_ini - 0.5, 'y': y}, 1: {'x': self.ch_end + 0.5, 'y': y}}

	def CreateLinesTLine(self):
		for lineh in self.horizontal_lines_diamond:
			self.horizontal_lines_diamond_tline.append(ro.TLine(lineh[0]['x'], lineh[0]['y'], lineh[1]['x'], lineh[1]['y']))
			self.horizontal_lines_diamond_tline[-1].SetLineColor(ro.kRed)
		for linev in self.vertical_lines_diamond:
			self.vertical_lines_diamond_tline.append(ro.TLine(linev[0]['x'], linev[0]['y'], linev[1]['x'], linev[1]['y']))
			self.vertical_lines_diamond_tline[-1].SetLineColor(ro.kRed)

	def ResetLines(self):
		self.horizontal_lines_diamond = []
		self.horizontal_lines_diamond_tline = []
		self.vertical_lines_diamond = []
		self.vertical_lines_diamond_tline = []

	def DrawLinesDiamond(self, name):
		self.DrawLines(name, type='diamond')

	def DrawLines(self, name, type):
		# ro.gStyle.SetOptStat('en')
		self.canvas[name].cd()
		if type == 'diamond':
			for lineh in self.horizontal_lines_diamond_tline:
				lineh.Draw('same')
			for linev in self.vertical_lines_diamond_tline:
				linev.Draw('same')

	def CreateTCutGs(self):
		self.CreateTCutGsDiamond()
		self.CreateTCutGsDiamondCenter()
		self.CreateGridText()

	def CreateTCutGsDiamond(self):
		def GetNumpyArraysX(coli):
			x0 = self.ch_ini - 0.5 + coli
			x1 = self.ch_ini + 0.5 + coli
			return np.array((x0, x0, x1, x1, x0), 'f8')
		def GetNumpyArraysY(rowi):
			y0 = self.row_info_diamond['0'] + rowi * self.row_info_diamond['pitch']
			y1 = self.row_info_diamond['0'] + (rowi + 1) * self.row_info_diamond['pitch']
			return np.array((y0, y1, y1, y0, y0), 'f8')
		for col in xrange(self.num_cols):
			self.tcutgs_diamond[col] = {}
			for row in xrange(self.row_info_diamond['num']):
				tempx = GetNumpyArraysX(col)
				tempy = GetNumpyArraysY(row)
				self.tcutgs_diamond[col][row] = ro.TCutG('cutg_dia_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
				self.tcutgs_diamond[col][row].SetNameTitle('cutg_dia_{c}_{r}'.format(c=col, r=row), 'cutg_dia_{c}_{r}'.format(c=col, r=row))
				self.tcutgs_diamond[col][row].SetVarX('diaChXPred')
				self.tcutgs_diamond[col][row].SetVarY('diaChYPred')
				self.tcutgs_diamond[col][row].SetLineColor(ro.kBlack)
				if col == self.num_cols -1 and row == self.row_info_diamond['num'] - 1:
					self.row_info_diamond['up'] = tempy[2]

	def CreateTCutGsDiamondCenter(self):
		def GetNumpyArraysX(coli):
			x0 = self.ch_ini - self.length_central_region/(2.0*self.col_pitch) + coli
			x1 = self.ch_ini + self.length_central_region/(2.0*self.col_pitch) + coli
			return np.array((x0, x0, x1, x1, x0), 'f8')

		def GetNumpyArraysY(rowi):
			y0 = self.row_info_diamond['0'] + rowi * self.row_info_diamond['pitch'] + self.row_info_diamond['pitch']/2.0 - self.length_central_region/2.0
			y1 = self.row_info_diamond['0'] + (rowi + 1) * self.row_info_diamond['pitch'] - self.row_info_diamond['pitch']/2.0 + self.length_central_region/2.0
			return np.array((y0, y1, y1, y0, y0), 'f8')

		x0i = self.col_pitch / 2.0 - self.length_central_region / 2.0
		x1i = self.col_pitch / 2.0 + self.length_central_region / 2.0
		y0i = self.row_info_diamond['pitch'] / 2.0 - self.length_central_region / 2.0
		y1i = self.row_info_diamond['pitch'] / 2.0 + self.length_central_region / 2.0
		tempi = np.array((x0i, x0i, x1i, x1i, x0i), 'f8')
		tempj = np.array((y0i, y1i, y1i, y0i, y0i), 'f8')
		self.tcutg_diamond_center = ro.TCutG('cutg_dia_center', 5, tempi, tempj)
		self.tcutg_diamond_center.SetNameTitle('cutg_dia_center', 'cutg_dia_center')
		self.tcutg_diamond_center.SetVarX('diaChXPred')
		self.tcutg_diamond_center.SetVarY('diaChYPred')
		self.tcutg_diamond_center.SetLineColor(ro.kViolet)

		for col in xrange(self.num_cols):
			self.tcutgs_diamond_center[col] = {}
			for row in xrange(self.row_info_diamond['num']):
				tempx = GetNumpyArraysX(col)
				tempy = GetNumpyArraysY(row)
				self.tcutgs_diamond_center[col][row] = ro.TCutG('cutg_dia_center_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
				self.tcutgs_diamond_center[col][row].SetNameTitle('cutg_dia_center_{c}_{r}'.format(c=col, r=row), 'cutg_center_dia_{c}_{r}'.format(c=col, r=row))
				self.tcutgs_diamond_center[col][row].SetVarX('diaChXPred')
				self.tcutgs_diamond_center[col][row].SetVarY('diaChYPred')
				self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)

	def CreateGridText(self):
		self.gridTextDiamond = ro.TH2F('gridText_diamond', 'gridText_diamond', int(np.floor(128.0 * self.bins_per_ch_x + 0.5) + 2), -0.5 - 1.0 / self.bins_per_ch_x, 127.5 + 1.0 / self.bins_per_ch_x, int(np.floor(256 * self.bins_per_ch_y + 0.5) + 2), self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'] - (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y), self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'] + (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y))
		x0, x1, y0, y1 = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		for col in xrange(0, self.num_cols):
			self.tcutgs_diamond[col][0].GetPoint(0, x0, y0)
			self.tcutgs_diamond[col][0].GetPoint(3, x1, y0)
			self.gridTextDiamond.Fill(np.mean((x0, x1)), y0[0]-0.1, (col + 0.01))
		for row in xrange(0, self.row_info_diamond['num']):
			self.tcutgs_diamond[0][row].GetPoint(0, x0, y0)
			self.tcutgs_diamond[0][row].GetPoint(1, x0, y1)
			self.gridTextDiamond.Fill(x0[0]-0.1, np.mean((y0, y1)), (row + 0.01))
		self.gridTextDiamond.SetMarkerSize(0.8)

	def DrawProfile2DDiamond(self, name, varz='clusterChargeN', cuts='', draw_top_borders=False, transp_ev=True):
		xmin, xmax, deltax, xname = -0.5, 127.5, 1.0/self.bins_per_ch_x, 'dia X ch'
		ymin, ymax, deltay, yname = self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'], float(self.row_info_diamond['pitch'])/self.bins_per_ch_y, 'sil pred Y [#mum]'
		if draw_top_borders:
			self.DrawProfile2D(name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, 'diaChXPred', 'diaChYPred', varz, 'PH[ADC]', cuts, transp_ev)
		else:
			self.DrawProfile2DNoTopBottomBorders(name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, 'diaChXPred', 'diaChYPred', varz, 'PH[ADC]', cuts, transp_ev)

	def DrawProfile2DNoTopBottomBorders(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, varz='clusterChargeN', zname='PH[ADC]', cuts='', transp_ev=True):
		list_cuts =['(({l}<=diaChYPred)&&(diaChYPred<={h}))'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['up'])]
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawProfile2D(name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, varz, zname, temp_cuts, transp_ev)

	def DrawProfile2D(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, varz='clusterChargeN', zname='PH[ADC]', cuts='', transp_ev=True):
		# ro.gStyle.SetOptStat('en')
		ro.TFormula.SetMaxima(100000)
		self.profile[name] = ro.TProfile2D('h_' + name, 'h_' + name, int(np.floor((xmax - xmin)/deltax + 0.5) + 2), xmin - deltax, xmax + deltax, int(np.floor((ymax - ymin)/deltay + 0.5) + 2), ymin - deltay, ymax + deltay)
		self.profile[name].GetXaxis().SetTitle(xname)
		self.profile[name].GetYaxis().SetTitle(yname)
		self.profile[name].GetZaxis().SetTitle(zname)
		self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		self.canvas[name].cd()
		list_cuts = ['transparentEvent'] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.trans_tree.Draw('{z}:{y}:{x}>>h_{n}'.format(z=varz, y=vary, x=varx, n=name), temp_cut, 'colz prof')
		ro.gPad.Update()
		SetDefault2DStats(self.profile[name])
		ro.TFormula.SetMaxima(1000)

	def DrawTCutGs(self, name, type):
		self.canvas[name].cd()
		# ro.gStyle.SetOptStat('en')
		ro.gStyle.SetPaintTextFormat(".0f")
		if type == 'diamond':
			self.gridTextDiamond.Draw('same TEXT0')
		if name in self.profile.keys():
			self.profile[name].Draw('same colz')
		elif name in self.histo.keys():
			self.histo[name].Draw('same colz')
		for col in xrange(0, self.num_cols):
			for row in xrange(0, self.row_info_diamond['num']):
				if type == 'diamond':
					self.tcutgs_diamond[col][row].Draw('same')
				elif type == 'centers':
					self.tcutgs_diamond_center[col][row].Draw('same')
		ro.gPad.Update()

	def GetOccupancyFromProfile(self, name):
		# ro.gStyle.SetOptStat('ne')
		name_occupancy = 'hit_map_' + name
		self.histo[name_occupancy] = self.profile[name].ProjectionXY('h_' + name_occupancy, 'B')
		self.histo[name_occupancy].SetTitle('h_' + name_occupancy)
		self.histo[name_occupancy].GetXaxis().SetTitle(self.profile[name].GetXaxis().GetTitle())
		self.histo[name_occupancy].GetYaxis().SetTitle(self.profile[name].GetYaxis().GetTitle())
		self.histo[name_occupancy].GetZaxis().SetTitle('entries')
		self.canvas[name_occupancy] = ro.TCanvas('c_' + name_occupancy, 'c_' + name_occupancy, 1)
		self.canvas[name_occupancy].cd()
		self.histo[name_occupancy].Draw('colz')
		ro.gPad.Update()
		SetDefault2DStats(self.histo[name_occupancy])

	def Draw2DHistoDiamond(self, name, cuts='', transp_ev=True):
		self.DrawHisto2D(name, -0.5, 127.5, 1.0 / (self.bins_per_ch_x), 'dia X ch', self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'],
						 float(self.row_info_diamond['pitch']) / (self.bins_per_ch_y), 'dia Y [#mum]', 'diaChXPred', 'diaChYPred', cuts, transp_ev)

	def DrawHisto2D(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, cuts='', transp_ev=True):
		# ro.gStyle.SetOptStat('en')
		ro.TFormula.SetMaxima(100000)
		self.histo[name] = ro.TH2F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5) + 2), xmin - deltax, xmax + deltax, int(np.floor((ymax - ymin) / deltay + 0.5) + 2), ymin - deltay, ymax + deltay)
		self.histo[name].GetXaxis().SetTitle(xname)
		self.histo[name].GetYaxis().SetTitle(yname)
		self.histo[name].GetZaxis().SetTitle('entries')
		self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		self.canvas[name].cd()
		list_cuts = ['transparentEvent'] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.trans_tree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), temp_cut, 'colz')
		ro.gPad.Update()
		SetDefault2DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def DrawPH(self, name, xmin, xmax, deltax, var='clusterChargeN', varname='PH[ADC]', cuts='', transp_ev=True, option='e'):
		ro.TFormula.SetMaxima(100000)
		# ro.gStyle.SetOptStat('neMmRruo')
		self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5)), xmin, xmax)
		self.histo[name].GetXaxis().SetTitle(varname)
		self.histo[name].GetYaxis().SetTitle('entries')
		list_cuts = ['transparentEvent'] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		if 'goff' not in option:
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		self.trans_tree.Draw('{v}>>h_{n}'.format(v=var, n=name), temp_cuts, option)
		if 'goff' not in option:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			SetDefault1DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def SelectGoodAndBadByThreshold(self, val=500, var='clusterChargeN'):
		for col in xrange(self.num_cols):
			for row in xrange(self.row_info_diamond['num']):
				# self.temph = ro.TH1F('temphrc', 'temphrc', 200, 0, 4000)
				self.trans_tree.Draw(var+'>>temphrc(200,0,4000)', 'transparentEvent&&({n})'.format(n=self.tcutgs_diamond[col][row].GetName()), 'goff')
				temph = ro.gDirectory.Get('temphrc')
				if temph.GetMean() > val:
					self.gridAreas.AddGoodAreas(col, row, self.tcutgs_diamond, self.tcutgs_diamond_center)
				else:
					self.gridAreas.AddBadAreas(col, row, self.tcutgs_diamond, self.tcutgs_diamond_center)
				temph.Reset('ICES')
				temph.Delete()
				del temph

	def DrawGoodAreasDiamond(self, name):
		self.DrawGoodAreas(name, type='diamond')

	def DrawGoodAreasDiamondCenters(self, name):
		self.DrawGoodAreas(name, type='centers')

	def DrawGoodAreas(self, name, type):
		# ro.gStyle.SetOptStat('en')
		self.canvas[name].cd()
		if type == 'diamond':
			for area in self.gridAreas.goodAreas_diamond:
				area.Draw('same')
		elif type == 'centers':
			for area in self.gridAreas.goodAreas_diamond_centers:
				area.Draw('same')

	def DrawBadAreasDiamond(self, name):
		self.DrawBadAreas(name, type='diamond')

	def DrawBadAreasDiamondCenters(self, name):
		self.DrawBadAreas(name, type='centers')

	def DrawBadAreas(self, name, type):
		# ro.gStyle.SetOptStat('en')
		self.canvas[name].cd()
		if type == 'diamond':
			for area in self.gridAreas.badAreas_diamond:
				area.Draw('same')
		elif type == 'centers':
			for area in self.gridAreas.badAreas_diamond_centers:
				area.Draw('same')

	def ResetHistos(self):
		for histo in self.histo.itervalues():
			histo.Delete()
		self.histo = {}

	def ResetProfiles(self):
		for profile in self.profile.itervalues():
			profile.Delete()
		self.profile = {}

	def ResetCanvas(self):
		for canvas in self.canvas.itervalues():
			canvas.Clear()
			canvas.Close()
		self.canvas = {}

	def ResetPlots(self):
		self.ResetHistos()
		self.ResetProfiles()
		self.ResetCanvas()

	def AddGoodAreas(self, col, row):
		self.gridAreas.AddGoodAreas(col, row, self.tcutgs_diamond, self.tcutgs_diamond_center)

	def AddBadAreas(self, col, row):
		self.gridAreas.AddBadAreas(col, row, self.tcutgs_diamond, self.tcutgs_diamond_center)

	def AddGoodAreasRow(self, row, coli=0, colf=0):
		self.gridAreas.AddGoodAreasRow(row, coli, colf, self.tcutgs_diamond, self.tcutgs_diamond_center)

	def AddGoodAreasCol(self, col, rowi=0, rowf=0):
		self.gridAreas.AddGoodAreasCol(col, rowi, rowf, self.tcutgs_diamond, self.tcutgs_diamond_center)

	def AddRemainingToBadAreas(self):
		self.gridAreas.AddRemainingToBadAreas(self.tcutgs_diamond, self.tcutgs_diamond_center)

	def RemoveFromGoodArea(self, col, row):
		self.gridAreas.RemoveFromGoodArea(col, row, self.tcutgs_diamond, self.tcutgs_diamond_center)

	def ResetAreas(self):
		self.gridAreas.ResetAreas()

	def DrawPHGoodAreas(self, name, var='clusterChargeN', cuts='', type='diamond', transp_ev=True):
		list_cuts = ['{n}'.format(n=self.gridAreas.goodAreasCutNames_diamond if type == 'diamond' else '')]
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cut, transp_ev)

	def DrawPHBadAreas(self, name, var='clusterChargeN', cuts='', type='diamond', transp_ev=True):
		list_cuts = ['{n}'.format(n=self.gridAreas.badAreasCutNames_diamond if type == 'diamond' else '')]
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cut, transp_ev)

	def DrawPHCentralRegion(self, name, var='clusterChargeN', cells='good', cuts='', transp_ev=True):
		list_cuts = ['{n}'.format(n=self.gridAreas.goodAreasCutNames_diamond_centers) if cells == 'good' else '{n}'.format(n=self.gridAreas.badAreasCutNames_diamond_centers) if cells == 'bad' else '({n}||{m})'.format(n=self.gridAreas.goodAreasCutNames_diamond_centers, m=self.gridAreas.badAreasCutNames_diamond_centers)]
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cuts, transp_ev)

	def DrawProfile2DDiamondChannelOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True):
		list_cuts = ['{n}'.format(n=self.gridAreas.goodAreasCutNames_diamond) if cells == 'good' else '{n}'.format(n=self.gridAreas.badAreasCutNames_diamond) if cells == 'bad' else '(1)']
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		rowpitch, y0, xoff = self.row_info_diamond['pitch'], self.row_info_diamond['0'], self.row_info_diamond['x_off']
		self.DrawProfile2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', y0 - np.floor(y0 / rowpitch + 0.5) * rowpitch, y0 + (256 - np.floor(y0 / rowpitch + 0.5)) * rowpitch,
						   float(rowpitch)/self.bins_per_ch_y, 'dia Y [#mum]', '((diaChXPred-{o})*{p})%{p}'.format(o=xoff, p=self.col_pitch), 'diaChYPred', var, 'PH[ADC]', temp_cuts, transp_ev)

	def DrawProfile2DDiamondRowOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True):
		y0, rowpitch, numrows, yoff = self.row_info_diamond['0'], self.row_info_diamond['pitch'], self.row_info_diamond['num'], self.row_info_diamond['y_off']
		list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=y0, h=y0 + rowpitch * numrows)]
		if cells == 'good':
			list_cuts.append(self.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts.append(self.gridAreas.badAreasCutNames_diamond)
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawProfile2D(name, -0.5, 127.5, self.cell_resolution, 'dia X ch', 0, rowpitch, self.cell_resolution, 'dia Y [#mum]', 'diaChXPred', '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=yoff, srp=int(10000*rowpitch)), var, 'PH[ADC]', temp_cuts, transp_ev)

	def DrawProfile2DDiamondCellOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True):
		y0, rowpitch, numrows, xoff, yoff = self.row_info_diamond['0'], self.row_info_diamond['pitch'], self.row_info_diamond['num'], self.row_info_diamond['x_off'], self.row_info_diamond['y_off']
		list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=y0, h=y0 + rowpitch * numrows)]
		if cells == 'good':
			list_cuts.append(self.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts.append(self.gridAreas.badAreasCutNames_diamond)
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawProfile2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', 0, rowpitch, self.cell_resolution, 'dia Y [#mum]', '((diaChXPred-{ox})*{p})%{p}'.format(ox=xoff, p=self.col_pitch), '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=yoff, srp=int(10000*rowpitch)), var, 'PH[ADC]', temp_cuts, transp_ev)

	def DrawHisto2DDiamondChannelOverlay(self, name, cells='all', cuts='', transp_ev=True):
		rowpitch, y0, xoff = self.row_info_diamond['pitch'], self.row_info_diamond['0'], self.row_info_diamond['x_off']
		list_cuts = []
		if cells == 'good':
			list_cuts.append(self.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts.append(self.gridAreas.badAreasCutNames_diamond)
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawHisto2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', y0 - np.floor(y0 / rowpitch + 0.5) * rowpitch, y0 + (256 - np.floor(y0 / rowpitch + 0.5)) * rowpitch,
						 float(rowpitch) / self.bins_per_ch_y, 'dia Y [#mum]', '((diaChXPred-{o})*{p})%{p}'.format(o=xoff, p=self.col_pitch), 'diaChYPred', temp_cuts, transp_ev)

	def DrawHisto2DDiamondRowOverlay(self, name, cells='all', cuts='', transp_ev=True):
		y0, rowpitch, numrows, yoff = self.row_info_diamond['0'], self.row_info_diamond['pitch'], self.row_info_diamond['num'], self.row_info_diamond['y_off']
		list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=y0, h=y0 + rowpitch * numrows)]
		if cells == 'good':
			list_cuts.append(self.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts.append(self.gridAreas.badAreasCutNames_diamond)
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawHisto2D(name, -0.5, 127.5, 1.0 / self.bins_per_ch_x, 'dia X ch', 0, rowpitch, self.cell_resolution, 'dia Y [#mum]', 'diaChXPred', '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=yoff, srp=int(10000*rowpitch)), temp_cuts, transp_ev)

	def DrawHisto2DDiamondCellOverlay(self, name, cells='all', cuts='', transp_ev=True):
		y0, rowpitch, numrows, xoff, yoff = self.row_info_diamond['0'], self.row_info_diamond['pitch'], self.row_info_diamond['num'], self.row_info_diamond['x_off'], self.row_info_diamond['y_off']
		list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=y0, h=y0 + rowpitch * numrows)]
		if cells == 'good':
			list_cuts.append(self.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts.append(self.gridAreas.badAreasCutNames_diamond)
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		self.DrawHisto2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', 0, rowpitch, self.cell_resolution, 'dia Y [#mum]', '((diaChXPred-{ox})*{p})%{p}'.format(ox=xoff, p=self.col_pitch), '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=yoff, srp=int(10000*rowpitch)), temp_cuts, transp_ev)

	def DrawTCutCentersInCellOverlay(self, name):
		self.canvas[name].cd()
		self.tcutg_diamond_center.Draw('same')

	def FitLanGaus(self, name, conv_steps=100, color=ro.kRed):
		self.canvas[name].cd()
		self.langaus[name] = LanGaus(self.histo[name])
		self.langaus[name].LanGausFit(conv_steps)
		lowbin, highbin = self.histo[name].FindFirstBinAbove(0, 1), self.histo[name].FindLastBinAbove(0, 1)
		xlow, xhigh = self.histo[name].GetBinLowEdge(lowbin), self.histo[name].GetBinLowEdge(highbin + 1)
		self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
		self.line[name].SetLineColor(ro.kViolet + 1)
		self.line[name].SetLineWidth(4)
		fitmean = self.langaus[name].fit.Mean(xlow, xhigh)
		self.langaus[name].fit.Draw('same')
		self.langaus[name].fit.SetLineColor(color)
		self.line[name].Draw('same')
		ro.gPad.Update()
		self.histo[name].FindObject('stats').SetOptFit(1)
		self.histo[name].FindObject('stats').SetX1NDC(0.6)
		self.histo[name].FindObject('stats').SetX2NDC(0.9)
		self.histo[name].FindObject('stats').SetY1NDC(0.6)
		self.histo[name].FindObject('stats').SetY2NDC(0.9)
		AddLineToStats(self.canvas[name], 'Mean_{Fit}', fitmean)
		self.histo[name].SetStats(0)
		self.canvas[name].Modified()
		ro.gPad.Update()
		print '{n}: <PH> = {f}'.format(n=name, f=fitmean)

	def DrawDoubleLangaus(self, name, name1, name2, color=ro.kBlack):
		if self.langaus.has_key(name1) and self.langaus.has_key(name2):
			langaus1 = self.langaus[name1].fit
			langaus2 = self.langaus[name2].fit
			self.doubleLangaus[name] = ro.TF1(name, self.TwoLanGaus, 0, 4000, 8)
			params1 = np.zeros(4, 'f8')
			params2 = np.zeros(4, 'f8')
			langaus1.GetParameters(params1)
			langaus2.GetParameters(params2)
			params12 = np.concatenate((params1, params2))
			self.doubleLangaus[name].SetNpx(1000)
			self.doubleLangaus[name].SetParameters(params12)
			self.doubleLangaus[name].SetParNames('Width1', 'MP1', 'Area1', 'GSigma1', 'Width2', 'MP2', 'Area2', 'GSigma2')
			lowbin, highbin = self.histo[name].FindFirstBinAbove(0, 1), self.histo[name].FindLastBinAbove(0, 1)
			xlow, xhigh = self.histo[name].GetBinLowEdge(lowbin), self.histo[name].GetBinLowEdge(highbin + 1)
			fitmean = self.doubleLangaus[name].Mean(xlow, xhigh)
			self.canvas[name].cd()
			self.doubleLangaus[name].Draw('same')
			self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
			self.line[name].SetLineColor(ro.kViolet + 1)
			self.line[name].SetLineWidth(4)
			self.line[name].Draw('same')
			self.doubleLangaus[name].SetLineColor(color)
			AddLineToStats(self.canvas[name], 'Mean_{Fit}', fitmean)
			self.histo[name].SetStats(0)
			self.canvas[name].Modified()
			ro.gPad.Update()
			print '{n}: <PH> = {f}'.format(n=name, f=fitmean)

	def TwoLanGaus(self, x, params):
		mpc1 = params[1] - self.mpshift * params[0]
		mpc2 = params[5] - self.mpshift * params[4]
		xlow1, xhigh1 = [x[0] + self.sigma_conv * i * params[3] for i in [-1, 1]]
		xlow2, xhigh2 = [x[0] + self.sigma_conv * i * params[7] for i in [-1, 1]]
		step1 = (xhigh1 - xlow1) / self.conv_steps
		step2 = (xhigh2 - xlow2) / self.conv_steps
		sums1 = 0
		sums2 = 0
		for i in xrange(1, int(np.ceil(self.conv_steps / 2.0 + 1))):
			xx1 = xlow1 + (i - 0.5) * step1
			xx2 = xlow2 + (i - 0.5) * step2
			fland1 = ro.TMath.Landau(xx1, mpc1, params[0]) / params[0]
			fland2 = ro.TMath.Landau(xx2, mpc2, params[4]) / params[4]
			sums1 += fland1 * ro.TMath.Gaus(x[0], xx1, params[3])
			sums2 += fland2 * ro.TMath.Gaus(x[0], xx2, params[7])
			xx1 = xhigh1 - (i - 0.5) * step1
			xx2 = xhigh2 - (i - 0.5) * step2
			fland1 = ro.TMath.Landau(xx1, mpc1, params[0]) / params[0]
			fland2 = ro.TMath.Landau(xx2, mpc2, params[4]) / params[4]
			sums1 += fland1 * ro.TMath.Gaus(x[0], xx1, params[3])
			sums2 += fland2 * ro.TMath.Gaus(x[0], xx2, params[7])
		return params[2] * step1 * sums1 / (np.sqrt(2 * np.pi, dtype='f8') * params[3]) + params[6] * step2 * sums2 / (np.sqrt(2 * np.pi, dtype='f8') * params[7])

	def SaveCanvasInlist(self, list):
		if not os.path.isdir('{d}/{r}/{sd}'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/{sd}'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir))
		for canvas in list:
			if self.canvas.has_key(canvas):
				self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.png'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))
				self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.root'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
	parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
	parser.add_option('-c', '--cellsize', dest='cellsize', type='int', default=50, help='cell size of the square 3D device')
	parser.add_option('-n', '--numstrips', dest='numstrips', type='int', default=2, help='Number of strips to use')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	# testnum = int(options.testnumber)
	numstrips = int(options.numstrips)
	cellsize = int(options.cellsize)

	tg = TransparentGrid(dir=dir, run=run, cellsize=cellsize)