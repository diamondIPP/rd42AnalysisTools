#!/usr/bin/env python
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import scipy.stats as sps
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
from copy import deepcopy
from glob import glob
from collections import OrderedDict
import os, sys, shutil
from Utils import *
import cPickle as pickle
from Langaus import LanGaus
from DiamondColumns import DiamondColumns
from GridAreas import GridAreas
from CutManager import CutManager
from PedestalCalculations import PedestalCalculations
from array import array
import time

# ph_bins_options = np.array((1, 2, 4, 5, 10, 16, 20, 25, 32, 40, 50, 80, 100, 125, 160, 200, 250, 400, 500, 800, 1000, 2000), 'uint16')
# ph_bins_options = np.array((32, 40, 50, 80, 100, 125, 160, 200, 250, 400, 500, 800, 1000, 2000, 4000), 'uint16')
# ph_bins_options = np.array((32, 40, 50, 80, 100, 125, 160, 200, 250, 400, 500), 'uint16')
ph_bins_options = np.array((25, 30, 32, 40, 50, 60, 64, 75, 80, 100, 120, 150, 160, 200, 240, 300, 320, 400, 480, 600), 'uint16')
eff_delta_options = np.array((8.2, 10, 10.25, 12.5, 16.4, 20, 20.5, 25, 32.8, 41, 50, 51.25, 82), 'float64')

class TransparentGrid:
	def __init__(self, dir='', run=25209, col_pitch=50.0):
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(99)
		ro.gStyle.SetOptStat(11)
		ro.gStyle.SetOptFit(111)
		ro.gStyle.SetPaintTextFormat(".0f")
		# ro.TFormula.SetMaxima(100000)
		self.num_parallel = 3
		self.run = run
		self.dir = os.path.abspath(os.path.expanduser(os.path.expandvars(dir)))
		self.trans_file = None
		self.trans_tree = None
		self.align_file = None
		self.align_obj = None
		self.pkl = None
		self.pkl_sbdir = 'default'
		self.dut = 'dia'
		self.loaded_pickle = False
		self.align_info = {'xoff': float(0), 'phi': float(0)}
		self.col_overlay_var = ''
		self.row_overlay_var = ''
		self.col_overlay_var2 = ''
		self.row_overlay_var2 = ''
		self.hit_factor = 0
		self.seed_factor = 0
		self.cluster_size, self.num_strips = 0, 0
		self.num_cols = 19 if col_pitch == 50 else 13
		self.ch_ini = 0
		self.ch_end = 127
		self.num_sides = 4
		self.threshold_criteria_in_sigmas = 2
		self.threshold = 0
		self.threshold_snr = 0
		self.phbins = 240
		self.phmin = 0
		self.phmax = 4800
		self.phbins_neg = 240
		self.phmin_neg = -2000
		self.phmax_neg = 2000
		self.neg_cut_snr = 410
		self.neg_cut_adc = 4100
		self.evs_before_sat_cut = 0
		self.evs_after_sat_cut = 0
		self.col_pitch = col_pitch
		self.cell_resolution = 0
		self.delta_offset_threshold = 0.01  # in mum
		self.saturated_ADC = 4095
		self.bias = 0
		self.row_cell_info_diamond = {'height': self.col_pitch, 'width': self.col_pitch, 'num_even': 0, 'num_odd': 0, '0_even': 0, '0_odd': 0, 'up_even': 12750, 'up_odd': 12750, 'x_off': 0, 'y_off': 0}
		self.bins_per_ch_x = int(1 + 2 * np.ceil(np.divide(self.col_pitch, 50.)).astype('int16'))
		self.bins_per_ch_y = int(1 + 2 * np.ceil(np.divide(self.col_pitch, 50.)).astype('int16'))
		self.length_central_region = 30 if self.col_pitch == 50 else 40 if self.col_pitch == 100 else 50
		self.conv_steps = 1000
		self.sigma_conv = 5
		self.efficiency_subdiv = 1
		self.xoffset, self.yoffset = 0, 0
		self.vertical_lines_diamond = []
		self.vertical_lines_diamond_tline = []
		self.horizontal_lines_diamond = []
		self.horizontal_lines_diamond_tline = []
		self.mpshift = -0.22278298
		self.noise_friend_buffer = 0
		self.canvas = {}
		self.line = {}
		self.profile = {}
		self.histo = {}
		self.graph = {}
		self.fits = {}
		self.names = []
		self.dia_cols = None
		self.tcutgs_diamond = {}
		self.tcutgs_diamond_center = {}
		self.gridAreas = None
		self.cuts_man = None
		self.temph = None
		self.langaus = {}
		self.gridTextDiamond = None
		self.doubleLangaus = {}
		self.lg1, self.lg2 = ro.TF1(), ro.TF1()
		self.CheckFoldersAndFiles()
		self.OpenFileAndGetTree()
		self.FindDiamondChannelLimits()
		self.list_neg_cuts_clusters = {}
		self.list_neg_cuts_noise = {}
		self.mean_ph_cell_dic = {'snr': {}, 'adc': {}}
		self.suffix = {'all': 'all', 'good': 'selected', 'bad': 'not_selected'}

		self.noise_varz = {'adc': 'diaChSignal', 'snr': 'diaChSignal/diaChPedSigmaCmc'}
		self.noise_friend_varz = {'adc': 'pedTree.diaChSignal', 'snr': 'pedTree.diaChSignal/pedTree.diaChPedSigmaCmc'}
		self.ph_adc_h_varz = {}
		self.ph_adc_ch_varz = {}
		self.ph_snr_h_varz = {}
		self.ph_snr_ch_varz = {}

		self.phN_adc_h_varz = {}
		self.phN_adc_ch_varz = {}
		self.phN_snr_h_varz = {}
		self.phN_snr_ch_varz = {}

		self.minz = {}
		self.maxz = {}

		self.ped_file = None
		self.ped_tree = None
		self.trash = []
		self.LoadPickle()

	def CheckFoldersAndFiles(self):
		if not os.path.isdir(self.dir):
			ExitMessage('The path to the directory "{d}" does not exist. Exiting...'.format(d=self.dir), code=os.EX_NOINPUT)
		if not os.path.isdir('{d}/{r}'.format(d=self.dir, r=self.run)):
			ExitMessage('There is no run {r} in the directory "{d}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)
		if not os.path.isfile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run)):
			ExitMessage('There is no transparent root file "transparent.{r}.root" in the directory "{d}/{r}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)

	def OpenFileAndGetTree(self, mode='READ'):
		mode = mode if mode in ['READ', 'UPDATE'] else 'READ'
		self.trans_file = ro.TFile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run), mode)
		self.trans_tree = self.trans_file.Get('transparentTree')
		if self.trans_tree.FindLeaf('numStrips'):
			self.num_strips = self.trans_tree.GetMaximum('numStrips')
			self.num_strips = RoundInt(self.num_strips)
		if self.trans_tree.FindLeaf('clusterSize'):
			self.cluster_size = self.trans_tree.GetMaximum('clusterSize')
			self.cluster_size = RoundInt(self.cluster_size)

	def OpenPedestalFileAndTree(self):
		if os.path.isfile('{d}/{r}/pedestalData.{r}.root'.format(d=self.dir, r=self.run)):
			self.ped_file = ro.TFile('{d}/{r}/pedestalData.{r}.root'.format(d=self.dir, r=self.run), 'READ')
			if self.ped_file.FindKey('pedestalTree'):
				self.ped_tree = self.ped_file.Get('pedestalTree')
				self.ped_tree.SetBranchStatus('*', 0)
				for branch in ['eventNumber', 'DiaADC']:
					self.ped_tree.SetBranchStatus(branch, 1)

	def FindDiamondChannelLimits(self):
		"""
		Automatically finds the upper and lower channels. The difference between the two limits will define the number of columns in the transparent grid
		:return:
		"""
		temph = ro.TH1F('temph', 'temph', 128, -0.5, 127.5)
		self.trans_tree.Draw('diaChXPred>>temph', 'transparentEvent', 'goff')
		self.ch_ini = int(temph.GetBinCenter(temph.FindFirstBinAbove(1, 1)))
		self.ch_end = int(temph.GetBinCenter(temph.FindLastBinAbove(1, 1)))
		self.num_cols = self.ch_end - self.ch_ini + 1

	def SavePickle(self, saveDefault=False):
		"""
		Save the pickle object inside the sub-directory pkl_sbdir
		:param saveDefault: If saveDefault is true, it will also create a subfolder for default
		:return:
		"""
		object_dic = {}
		object_dic['row_cell_info_diamond'] = self.row_cell_info_diamond
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
		object_dic['neg_cut_snr'] = self.neg_cut_snr
		object_dic['neg_cut_adc'] = self.neg_cut_adc
		object_dic['col_pitch'] = self.col_pitch
		object_dic['cell_resolution'] = self.cell_resolution
		object_dic['bins_per_ch_x'] = self.bins_per_ch_x
		object_dic['bins_per_ch_y'] = self.bins_per_ch_y
		object_dic['length_central_region'] = self.length_central_region
		object_dic['conv_steps'] = self.conv_steps
		object_dic['sigma_conv'] = self.sigma_conv
		object_dic['mpshift'] = self.mpshift
		object_dic['efficiency_subdiv'] = self.efficiency_subdiv
		object_dic['saturated_ADC'] = self.saturated_ADC
		object_dic['bias'] = self.bias
		object_dic['events_after_saturation_cut'] = self.evs_after_sat_cut
		object_dic['events_before_saturation_cut'] = self.evs_before_sat_cut
		object_dic['num_parallel'] = self.num_parallel
		object_dic['hit_factor'] = self.hit_factor
		object_dic['seed_factor'] = self.seed_factor
		object_dic['num_strips'] = self.num_strips
		object_dic['num_sides'] = self.num_sides
		object_dic['dut'] = self.dut

		if not saveDefault:
			if not os.path.isdir('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
				os.makedirs('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
			picklepath = '{d}/{r}/{s}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)

			pickle.dump(object_dic, open(picklepath, 'wb'))

		else:
			if not os.path.isdir('{d}/{r}/default'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
				os.makedirs('{d}/{r}/default'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
				picklepath = '{d}/{r}/default/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
				pickle.dump(object_dic, open(picklepath, 'wb'))

		print 'Saved pickle :D'

	def LoadPickle(self):
		"""
		This method tries to read a pickle file in the subdirectory self.pkl_sbdir with configuration data for setting up the transparent grid
		If it is successful, self.loaded_pickle is True, otherwise False
		:return:
		"""
		picklepath = '{d}/{r}/{s}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
		if os.path.isfile(picklepath):
			with open(picklepath, 'rb') as pkl:
				self.pkl = pickle.load(pkl)
				self.UnfoldPickle()
				self.loaded_pickle = True
		else:
			print 'Coult not find the pickle: {d}/{r}/{s}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)
			self.loaded_pickle = False

	def UnfoldPickle(self):
		"""
		This Method takes the loaded pickle object and unfolds it into the transparent grid variables
		:return:
		"""
		if 'row_cell_info_diamond' in self.pkl.keys():
			keys = self.row_cell_info_diamond.keys()
			for key in keys:
				if key in self.pkl['row_cell_info_diamond'].keys():
					self.row_cell_info_diamond[key] = self.pkl['row_cell_info_diamond'][key]

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
		if 'neg_cut_snr' in self.pkl.keys():
			self.neg_cut_snr = self.pkl['neg_cut_snr']
		if 'neg_cut_adc' in self.pkl.keys():
			self.neg_cut_adc = self.pkl['neg_cut_adc']
		if 'col_pitch' in self.pkl.keys():
			self.col_pitch = self.pkl['col_pitch']
		if 'cell_resolution' in self.pkl.keys():
			self.cell_resolution = self.pkl['cell_resolution']
		if 'saturated_ADC' in self.pkl.keys():
			self.saturated_ADC = self.pkl['saturated_ADC']
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
		if 'efficiency_subdiv' in self.pkl.keys():
			self.efficiency_subdiv = self.pkl['efficiency_subdiv']
		if 'bias' in self.pkl.keys():
			self.bias = self.pkl['bias']
		if 'events_before_saturation_cut' in self.pkl.keys():
			self.evs_before_sat_cut = self.pkl['events_before_saturation_cut']
		if 'events_after_saturation_cut' in self.pkl.keys():
			self.evs_after_sat_cut = self.pkl['events_after_saturation_cut']
		if 'num_parallel' in self.pkl.keys():
			self.num_parallel = self.pkl['num_parallel']
		if 'hit_factor' in self.pkl.keys():
			self.hit_factor = self.pkl['hit_factor']
		if 'seed_factor' in self.pkl.keys():
			self.seed_factor = self.pkl['seed_factor']
		if 'num_strips' in self.pkl.keys():
			self.num_strips = self.pkl['num_strips']
		if 'num_sides' in self.pkl.keys():
			self.num_sides = self.pkl['num_sides']
		if 'dut' in self.pkl.keys():
			self.dut = self.pkl['dut']

	def CreateGridAreas(self):
		"""
		This method creates GridAreas object which will contain the information of the selected areas, not selected areas in the transparent grid
		:return:
		"""
		if self.gridAreas:
			self.gridAreas = None
		self.gridAreas = GridAreas(self.num_cols, self.row_cell_info_diamond['num_even'], self.row_cell_info_diamond['num_odd'], self.run)
		self.SetOverlayVariables()

	def SetOverlayVariables(self):
		"""
		This Method sets the variables for plots which overlay columns, rows or cells. Each time there is a change in self.xoffset or self.yoffset, it must be called to apply the changes
		:return:
		"""
		col_fact = 10000  # ~0.1nm
		row_fact = 10000  # ~0.1nm
		self.col_overlay_var = '(((diaChXPred-{ox})*({p}*{cf}))%({p}*{cf}))/{cf}'.format(ox=self.row_cell_info_diamond['x_off'], p=self.col_pitch, cf=col_fact)
		self.row_overlay_var = '(((diaChYPred-{oy})*{rf})%({rp}*{rf}))/{rf}'.format(oy=self.row_cell_info_diamond['y_off'], rp=self.row_cell_info_diamond['height'], rf=row_fact)
		self.col_overlay_var2 = '{cp}*((x0-{xo})-TMath::Floor(0.5+(x0-{xo})))'.format(cp=self.row_cell_info_diamond['width'], xo=self.xoffset)
		self.row_overlay_var2 = '{rp}*((y0-{yo})/{rp}-TMath::Floor(0.5+(y0-{yo})/{rp}))'.format(rp=self.row_cell_info_diamond['height'], yo=self.yoffset)

	def FindPickleValues(self, do_offset_plots=True, do_binning=False):
		"""
		Automatically find some pickle parameters that could be close to the real values. Should only be used at the beginning of the analysis to save some time
		:param do_offset_plots: If true, it will show the offset plots used to calculate the xoffset and yoffset.
		:param do_binning: if true, it will automatically calculate the binning parameters for ph 1d plots and the cell resolution for cell overlay plots
		:return:
		"""
		self.FindUpperAndLowerLimitsInY()
		self.SavePickle()
		if do_binning:
			self.FindBinningAndResolution()
			self.SavePickle()
		if self.trans_tree.GetFriend('cell{x}X{y}Y'.format(x=self.row_cell_info_diamond['x_off'], y=self.row_cell_info_diamond['y_off'])):
			self.FindXandYOffests(do_plot=do_offset_plots)
		self.SavePickle()

	def FindUpperAndLowerLimitsInY(self, cut=''):
		"""
		This Method tries to find the lower and higher limits in the Y axis for the detector. t will find even and odd independently. If the tree still has no friend with cells information, it will use the predicted channel position from the transparent tree.
		If it only has one column, then the values of odd and even are equal. If it is made of rectangles with equal number of rows, then it will average the values for even and odd.
		:param cut: Is a string with the cut that could be applied in the Draw methods. If no cuts are applied (the default) it is equal to ''
		:return:
		"""
		# Get the profile 2D map of the clusterChargeN
		col_types = ['even', 'odd']
		for col_type in col_types:
			if col_type == 'even' or self.num_cols > 1:
				if self.trans_tree.GetFriend('cells{x}X{y}Y'.format(x=self.row_cell_info_diamond['x_off'], y=self.row_cell_info_diamond['y_off'])):
					cut_typ = self.cuts_man.even_cols if col_type == 'even' else self.cuts_man.odd_cols
				else:
					cut_typ = self.cuts_man.even_pred_chs if col_type == 'even' else self.cuts_man.odd_pred_chs
				cuts = self.cuts_man.ConcatenateCuts([cut, cut_typ])
				self.DrawProfile2DDiamond('vertical_limits_profile_' + col_type, 'clusterChargeN', cuts=cuts, transp_ev=True)
				xbinmin, xbinmax = int(self.profile['vertical_limits_profile_' + col_type].GetXaxis().FindBin(self.ch_ini - 0.5)), int(self.profile['vertical_limits_profile_' + col_type].GetXaxis().FindBin(self.ch_ini - 0.5) + self.num_cols * self.bins_per_ch_x - 1)
				self.canvas['vertical_limits_profile_{t}_py'.format(t=col_type)] = ro.TCanvas('c_vertical_limits_profile_{t}_py'.format(t=col_type), 'c_vertical_limits_profile_{t}_py'.format(t=col_type), 1)
				self.canvas['vertical_limits_profile_{t}_py'.format(t=col_type)].cd()
				self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)] = self.profile['vertical_limits_profile_' + col_type].ProjectionY('h_vertical_limits_profile_{t}_py'.format(t=col_type), xbinmin, xbinmax, 'e hist')
				minbiny, maxbiny = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindFirstBinAbove(), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindLastBinAbove()
				for biny in xrange(maxbiny, int(self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetNbins())):
					if self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(biny) != 0:
						maxbiny = biny
				miny, maxy = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinLowEdge(minbiny), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinLowEdge(maxbiny + 1)
				self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().SetRangeUser(miny, maxy)
				if self.row_cell_info_diamond['num_' + col_type] != 0:
					# The number of cells on this type of column has been already established
					func = ro.TF1('box_fcn_' + col_type, '[0]*(TMath::Erf((x-([3]-{p}))/[1])+1)/2-[2]*(TMath::Erf((x-[3])/[4])+1)/2+[5]'.format(p=self.row_cell_info_diamond['num_' + col_type] * self.row_cell_info_diamond['height']), miny, maxy)
					func.SetNpx(10000)
					zmin, zmax = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetMinimum(), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetMaximum()
					y1bin, y2bin = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindFirstBinAbove((zmin + zmax) / 2.0), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindLastBinAbove((zmin + zmax) / 2.0) + 1
					y1, y2 = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinCenter(y1bin), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinCenter(y2bin)
					z0, z1, z2 = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((minbiny))), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((y1bin + y2bin) / 2.0)), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((maxbiny)))
					func.SetParLimits(0, abs(z1 - z0) / 100.0, 10.0 * abs(z1 - z0))
					func.SetParLimits(1, 0.1, 50)
					func.SetParLimits(2, abs(z1 - z2) / 100.0, 10.0 * abs(z1 - z2))
					func.SetParLimits(3, y2 - 500, y2 + 500)
					func.SetParLimits(4, 0.1, 50)
					func.SetParLimits(5, -1000.0 * abs(z0), 1000.0 * abs(z0))
					params = np.array((abs(z1 - z0), 20, abs(z1 - z2), y2, 20, z0), 'float64')
					func.SetParameters(params)
					fit_prof_proj_y = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].Fit('box_fcn_' + col_type, 'QIEBMSL', 'goff', self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinLowEdge(int((minbiny))), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinLowEdge(int((maxbiny))))
					params = np.array((fit_prof_proj_y.Parameter(0), fit_prof_proj_y.Parameter(1), fit_prof_proj_y.Parameter(2), fit_prof_proj_y.Parameter(3), fit_prof_proj_y.Parameter(4), fit_prof_proj_y.Parameter(5)), 'float64')
					func.SetParameters(params)
					fit_prof_proj_y = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].Fit('box_fcn_' + col_type, 'QIEBMSL', 'goff', (miny), (maxy))
					self.row_cell_info_diamond['0_' + col_type] = fit_prof_proj_y.Parameter(3) - self.row_cell_info_diamond['height'] * self.row_cell_info_diamond['num_' + col_type]
					self.row_cell_info_diamond['up_' + col_type] = fit_prof_proj_y.Parameter(3)
				else:
					# the number of cells on this type of column has not been stablishe. It will be established by fitting
					func = ro.TF1('box_fcn_' + col_type, '[0]*(TMath::Erf((x-([3]-[6]*{p}))/[1])+1)/2-[2]*(TMath::Erf((x-[3])/[4])+1)/2+[5]'.format(p=self.row_cell_info_diamond['height']), miny, maxy)
					func.SetNpx(10000)
					zmin, zmax = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetMinimum(), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetMaximum()
					y1bin, y2bin = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindFirstBinAbove((zmin + zmax) / 2.0), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].FindLastBinAbove((zmin + zmax) / 2.0) + 1
					y1, y2 = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinCenter(y1bin), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetXaxis().GetBinCenter(y2bin)
					z0, z1, z2 = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((minbiny))), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((y1bin + y2bin) / 2.0)), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinContent(int((maxbiny)))
					func.SetParLimits(0, abs(z1 - z0) / 100.0, 10.0 * abs(z1 - z0))
					func.SetParLimits(1, 0.1, 50)
					func.SetParLimits(2, abs(z1 - z2) / 100.0, 10.0 * abs(z1 - z2))
					func.SetParLimits(3, y2 - 500, y2 + 500)
					func.SetParLimits(4, 0.1, 50)
					func.SetParLimits(5, -1000.0 * abs(z0), 1000. * abs(z0))
					func.SetParLimits(6, max(1, abs(y2 - y1) / self.row_cell_info_diamond['height'] - 5), abs(y2 - y1) / self.row_cell_info_diamond['height'] + 5)
					params = np.array((abs(z1 - z0), 20, abs(z1 - z2), y2, 20, z0, abs(y2 - y1) / self.row_cell_info_diamond['height']), 'float64')
					func.SetParameters(params)
					fit_prof_proj_y = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].Fit('box_fcn_' + col_type, 'QIEBMSL', 'goff', self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinLowEdge(int((minbiny))), self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].GetBinLowEdge(int((maxbiny))))
					params = np.array((fit_prof_proj_y.Parameter(0), fit_prof_proj_y.Parameter(1), fit_prof_proj_y.Parameter(2), fit_prof_proj_y.Parameter(3), fit_prof_proj_y.Parameter(4), fit_prof_proj_y.Parameter(5), fit_prof_proj_y.Parameter(6)), 'float64')
					func.SetParameters(params)
					fit_prof_proj_y = self.histo['vertical_limits_profile_{t}_py'.format(t=col_type)].Fit('box_fcn_' + col_type, 'QIEBMSL', 'goff', (miny), (maxy))
					self.row_cell_info_diamond['num_' + col_type] = np.floor(fit_prof_proj_y.Parameter(6))
					self.row_cell_info_diamond['up_' + col_type] = fit_prof_proj_y.Parameter(3)
					self.row_cell_info_diamond['0_' + col_type] = fit_prof_proj_y.Parameter(3) - self.row_cell_info_diamond['height'] * self.row_cell_info_diamond['num_' + col_type]

		if self.num_cols == 1:
			self.row_cell_info_diamond['0_odd'] = self.row_cell_info_diamond['0_even']
			self.row_cell_info_diamond['up_odd'] = self.row_cell_info_diamond['up_even']
			self.row_cell_info_diamond['num_odd'] = self.row_cell_info_diamond['num_even']

		if self.num_sides == 4 and self.row_cell_info_diamond['num_even'] == self.row_cell_info_diamond['num_odd']:
			# if the cells are rectangles and the number of cells for even and odd columns is the same, then it is the old type and the values for even and odd will be averaged
			y_0 = np.array([self.row_cell_info_diamond['0_even'], self.row_cell_info_diamond['0_odd']]).mean()
			y_up = np.array([self.row_cell_info_diamond['up_even'], self.row_cell_info_diamond['up_odd']]).mean()
			for col_type in col_types:
				self.row_cell_info_diamond['0_' + col_type] = y_0
				self.row_cell_info_diamond['up_' + col_type] = y_up

	def ShiftVerticalLimits(self, yoff=0, updatePickle=False):
		"""
		Used to quickly shift the vertical variables in the pickle by a fixed quantity 'yoff'
		:param yoff: value in mum to shift the vertical variables
		:param updatePickle: If true, the pickle will be saved after the modification
		:return:
		"""
		self.row_cell_info_diamond['0_even'] += yoff
		self.row_cell_info_diamond['0_odd'] += yoff
		self.row_cell_info_diamond['up_odd'] += yoff
		self.row_cell_info_diamond['up_even'] += yoff
		if updatePickle:
			self.SavePickle()

	def SaveLocalOffestInPickle(self):
		"""
		Used to quickly transfer the offset results from FindXOffset and FindYOffset into the pickle offset variables to be saved
		:return:
		"""
		self.row_cell_info_diamond['x_off'] += self.xoffset
		self.row_cell_info_diamond['y_off'] += self.yoffset
		self.SavePickle()

	def FindBinningAndResolution(self):
		"""
		Method used to automatically find a binning suitable for analysis. It sets the delta_adc (phbins) and the cell_resolution
		It should be avoided if a desired binning is desired. e.g. to compare multiple runs
		:return:
		"""
		if self.gridAreas:
			if len(self.gridAreas.goodAreas_index) == 0:
				self.ResetAreas()
				self.CreateGridAreas()
				self.CreateTCutGs()
				if self.threshold != 0:
					self.SelectGoodAndBadByThreshold()
				else:
					self.FindThresholdCutFromCells()
					self.SelectGoodAndBadByThreshold()
				self.AddRemainingToBadAreas()
				self.gridAreas.SimplifyGoodAndBadAreas()
				self.cuts_man.SetCells(selection=self.gridAreas.goodAreasCutNames_simplified_diamond, not_selection=self.gridAreas.badAreasCutNames_simplified_diamond)
		else:
			self.CreateGridAreas()
			self.CreateTCutGs()
			self.SelectGoodAndBadByThreshold()
			if self.threshold != 0:
				self.SelectGoodAndBadByThreshold()
			else:
				self.FindThresholdCutFromCells()
				self.SelectGoodAndBadByThreshold()
			self.AddRemainingToBadAreas()
			self.gridAreas.SimplifyGoodAndBadAreas()
			self.cuts_man.SetCells(selection=self.gridAreas.goodAreasCutNames_simplified_diamond, not_selection=self.gridAreas.badAreasCutNames_simplified_diamond)
		self.DrawPHGoodAreas('binning_temp', 'clusterCharge1')
		histo_entries = float(self.histo['binning_temp'].GetEntries())
		temp = np.abs(np.subtract(ph_bins_options, histo_entries * 240.0 / 4800.0, dtype='float32'), dtype='float32')
		self.phbins = ph_bins_options[temp.argmin()]
		self.phbins_neg = ph_bins_options[temp.argmin()]
		cell_bins = min(25, max(11, RoundInt(np.sqrt(histo_entries * ((self.col_pitch / 2.0) ** 2) / 4800.0, dtype='float64'))))
		self.cell_resolution = np.divide(self.col_pitch, cell_bins, dtype='float64') if cell_bins % 2 == 1 else np.divide(50.0, cell_bins + 1, dtype='float64')
		self.DrawPHGoodAreas('binning_temp', 'clusterCharge1')

	def FindXandYOffests(self, factor=0.2, do_plot=False, cells='all'):
		"""
		Method to find offsets in x and y
		:param factor: factor used for the increment in each iteration. It should be > 0 and it is recommended to be <= 1
		:param do_plot: Show plots
		:param cells: Choose which cells to use for the fine alignment
		:return:
		"""
		self.FindXOffset(factor, do_plot, cells)
		self.FindYOffset(factor, do_plot, cells)

	def ShiftHalfXandYOffsets(self, n=1):
		"""
		Method used to shift the cell in X and in Y by half a cell. It is usefull to check fine-alignment
		:param n: integer which states how many odd-integer halves of the cell must be shifted
		:return:
		"""
		tempn = RoundInt(n)
		tempn = tempn - 0.5 if tempn > 0 else tempn + 0.5 if tempn < 0 else 0
		if tempn == 0:
			print 'must give an integer different from 0'
			return
		# self.row_cell_info_diamond['x_off'] += tempn
		# self.row_cell_info_diamond['y_off'] += tempn * self.row_cell_info_diamond['height']
		self.xoffset += tempn * self.row_cell_info_diamond['width'] / float(self.col_pitch)
		self.yoffset += tempn * self.row_cell_info_diamond['height']
		self.SetOverlayVariables()

	def ShiftHalfXOffset(self, n=1):
		"""
		Method used to shift the cell in X by half a cell. It is usefull to check fine-alignment
		:param n: integer which states how many odd-integer halves of the cell must be shifted
		:return:
		"""
		tempn = RoundInt(n)
		tempn = tempn - 0.5 if tempn > 0 else tempn + 0.5 if tempn < 0 else 0
		if tempn == 0:
			print 'must give an integer different from 0'
			return
		# self.row_cell_info_diamond['x_off'] += tempn
		self.xoffset += tempn * self.row_cell_info_diamond['width'] / float(self.col_pitch)
		self.SetOverlayVariables()

	def ShiftHalfYOffset(self, n=1):
		"""
		Method used to shift the cell in Y by half a cell. It is usefull to check fine-alignment
		:param n: integer which states how many odd-integer halves of the cell must be shifted
		:return:
		"""
		tempn = RoundInt(n)
		tempn = tempn - 0.5 if tempn > 0 else tempn + 0.5 if tempn < 0 else 0
		if tempn == 0:
			print 'must give an integer different from 0'
			return
		# self.row_cell_info_diamond['y_off'] += tempn * self.row_cell_info_diamond['height']
		self.yoffset += tempn * self.row_cell_info_diamond['height']
		self.SetOverlayVariables()

	def FindXOffset(self, factor=0.2, do_plot=False, cells='all'):
		"""
		Find the offset on X
		:param factor: factor used for the increment in each iteration. It should be > 0 and it is recommended to be <= 1
		:param do_plot: Show plots while calculating
		:param cells: Choose which cells to use for the fine alignment. Could be 'good', 'bad' , 'all'
		:return:
		"""
		plot_option = 'prof colz' if do_plot else 'prof goff'
		self.ShiftHalfXOffset(2)
		delta_x = self.row_cell_info_diamond['width']
		# proj_width = self.row_cell_info_diamond['height'] / 5.0  # in mum
		proj_width = self.row_cell_info_diamond['height']  # in mum
		proj_bins = RoundInt(float(proj_width) / self.cell_resolution)
		proj_bins = proj_bins if proj_bins % 2 == 1 else proj_bins + 1
		proj_low = int(RoundInt(RoundInt(self.row_cell_info_diamond['height'] / float(self.cell_resolution), 'f8') / 2.0, 'f8') - (RoundInt(proj_bins / 2.0, 'f8') - 1) + 1)
		proj_high = proj_low + proj_bins - 1
		iteration = 0
		min_delta = self.row_cell_info_diamond['width']
		x_off_shifted, x_min = 0.0, 0.0
		while (abs(delta_x) > self.delta_offset_threshold and iteration < 100) or iteration < 10:
			self.xoffset -= np.divide(delta_x, self.row_cell_info_diamond['width'], dtype='f8') * (np.exp(-iteration) * (1 - factor) + factor)
			# self.row_cell_info_diamond['x_off'] -= np.divide(delta_x, self.col_pitch, dtype='float64') * (np.exp(-iteration) * (1 - factor) + factor)
			self.SetOverlayVariables()
			self.DrawProfile2DDiamondCellOverlay('x_off_alignment', 'clusterCharge1', cells, plot_option=plot_option)
			# h_proj_x = self.profile['x_off_alignment'].ProjectionX('x_off_alignment_px', proj_low, proj_high)
			h_proj_x = self.profile['x_off_alignment'].ProjectionX('x_off_alignment_px', proj_low, proj_high, 'e hist')
			h_proj_x.GetXaxis().SetRangeUser(-self.row_cell_info_diamond['width'] / 2.0 + 0.1, self.row_cell_info_diamond['width'] / 2.0 - 0.1)
			minx = h_proj_x.GetBinCenter(h_proj_x.GetMinimumBin())
			minx = minx if abs(minx - self.row_cell_info_diamond['width'] / 2.0) > self.cell_resolution * 1.5 else self.row_cell_info_diamond['width'] / 2.0
			bins_fit_range = max(RoundInt(2.0 * 7 / self.cell_resolution, 'f8'), 3)
			fit_px = h_proj_x.Fit('pol2', 'QEFSN', '', max(minx - bins_fit_range * self.cell_resolution / 2.0, -self.row_cell_info_diamond['width'] / 2.0), min(minx + bins_fit_range * self.cell_resolution / 2.0, self.row_cell_info_diamond['width'] / 2.0))
			x_min = -fit_px.Parameter(1) / (2 * fit_px.Parameter(2))
			delta_x = 0 - x_min
			if abs(delta_x) < min_delta:
				min_delta = abs(delta_x)
				x_off_shifted = self.xoffset
			iteration += 1
			print iteration, x_min, delta_x
		print 'final', min_delta
		# self.row_cell_info_diamond['x_off'] = x_off_shifted - 0.5
		self.xoffset = x_off_shifted - 0.5 * self.row_cell_info_diamond['width'] / float(self.col_pitch)
		self.SetOverlayVariables()

	def FindYOffset(self, factor=0.2, do_plot=False, cells='all'):
		"""
		Find the offset on Y
		:param factor: factor used for the increment in each iteration. It should be > 0 and it is recommended to be <= 1
		:param do_plot: Show plots while calculating
		:param cells: Choose which cells to use for the fine alignment. Could be 'good', 'bad' , 'all'
		:return:
		"""
		plot_option = 'prof colz' if do_plot else 'prof goff'
		self.ShiftHalfYOffset(2)
		delta_y = self.row_cell_info_diamond['height']
		proj_width = self.row_cell_info_diamond['width']  # in mum
		# proj_width = self.col_pitch / 5.0  # in mum
		proj_bins = RoundInt(float(proj_width) / self.cell_resolution)
		proj_bins = proj_bins if proj_bins % 2 == 1 else proj_bins + 1
		proj_low = int(RoundInt(RoundInt(float(self.col_pitch) / self.cell_resolution, 'f8') / 2.0, 'f8') - (RoundInt(proj_bins / 2.0, 'f8') - 1) + 1)
		proj_high = proj_low + proj_bins - 1
		iteration = 0
		min_delta = self.row_cell_info_diamond['height']
		y_off_shifted, y_min = 0.0, 0.0
		while (abs(delta_y) > self.delta_offset_threshold and iteration < 100) or iteration < 10:
			self.yoffset -= delta_y * (np.exp(-iteration) * (1 - factor) + factor)
			# self.row_cell_info_diamond['y_off'] -= delta_y * (np.exp(-iteration) * (1 - factor) + factor)
			self.SetOverlayVariables()
			self.DrawProfile2DDiamondCellOverlay('y_off_alignment', 'clusterCharge1', cells, plot_option=plot_option)
			h_proj_y = self.profile['y_off_alignment'].ProjectionY('y_off_alignment_py', proj_low, proj_high, 'e hist')
			h_proj_y.GetXaxis().SetRangeUser(-self.row_cell_info_diamond['height'] / 2.0 + 0.1, self.row_cell_info_diamond['height'] / 2.0 - 0.1)
			miny = h_proj_y.GetBinCenter(h_proj_y.GetMinimumBin())
			miny = miny if abs(miny - self.row_cell_info_diamond['height'] / 2.0) > self.cell_resolution * 1.5 else self.row_cell_info_diamond['height'] / 2.0
			fit_py = h_proj_y.Fit('pol2', 'QEFSN', '', max(miny - 2.5 * self.cell_resolution - 1e-12, -self.row_cell_info_diamond['height'] / 2.0), min(miny + 2.5 * self.cell_resolution + 1e-12, self.row_cell_info_diamond['height'] / 2.0))
			y_min = -fit_py.Parameter(1) / (2 * fit_py.Parameter(2))
			delta_y = 0 - y_min
			if abs(delta_y) < min_delta:
				min_delta = abs(delta_y)
				y_off_shifted = self.yoffset
			iteration += 1
			print iteration, y_min, delta_y
		print 'final', min_delta
		# self.row_cell_info_diamond['y_off'] = y_off_shifted - self.row_cell_info_diamond['height'] / 2.0
		self.yoffset = y_off_shifted - 0.5 * self.row_cell_info_diamond['height']
		self.SetOverlayVariables()

	def SetupCutManager(self):
		"""
		Deletes a previous CutManager object if it existed, and creates a new one with the current parameters
		:return:
		"""
		if self.cuts_man:
			self.cuts_man = None
		self.cuts_man = CutManager(self.trans_tree, self.num_strips, self.cluster_size, self.saturated_ADC, self.neg_cut_adc, self.neg_cut_snr)

	def CreateTCutGs(self):
		"""
		Creates the TCutGs for analysis. This is done by creating dia_cols variable which is an object of the class DiamondColumns.
		This will contain all the information of individual columns and cells in the transparent grid. It is used to be able to use different cell-geometries such as hexagons or rectangles
		:return:
		"""
		self.CreateTCutGsDiamond()
		self.CreateGridText()
		self.SetupCutManager()

	def CreateTCutGsDiamond(self):
		"""
		Creates the DiamondColumns object with the offsets in x and y specified by the pickle
		:return:
		"""
		self.CreateTCutGsDiaCols(self.row_cell_info_diamond['x_off'], self.row_cell_info_diamond['y_off'])

	def CreateTCutGsDiaCols(self, xoff=0, yoff=0):
		"""
		Creates the DiamondColumns object with the specified offsets entered as parameters.
		:param xoff: is a value which indicates how much of the column should be shifted horizontally and in what direction. i.e. a value of 0.5 will indicate a shifting of the center of the columns by half the width of the column
		:param yoff: is a value which indicates how many um should be shifted vertically and in what direction. i.e. a value of 50 will indicate a shifting of the center of the rows by 50um upwards
		:return:
		"""
		if self.dia_cols:
			self.dia_cols = None
		self.dia_cols = DiamondColumns(self.num_cols, self.row_cell_info_diamond['height'], self.num_sides, self.run, self.row_cell_info_diamond['width'] / float(self.col_pitch))
		for col in xrange(self.num_cols):
			if RoundInt(self.ch_ini + col + xoff) % 2 == 0:
				self.dia_cols.SetupColumns(col, self.row_cell_info_diamond['num_even'], self.ch_ini + col + xoff, self.row_cell_info_diamond['0_even'] + yoff)
				self.dia_cols.cols[col].SetCellsInColumn()
			else:
				self.dia_cols.SetupColumns(col, self.row_cell_info_diamond['num_odd'], self.ch_ini + col + xoff, self.row_cell_info_diamond['0_odd'] + yoff)
				self.dia_cols.cols[col].SetCellsInColumn()

	def CreateTCutGSymmetricRectangle(self, percentage=80):
		"""
		TODO
		:param percentage:
		:return:
		"""
		self.gridAreas.CreateRectangleSymmetricCentralRegion(percentage, self.col_pitch, self.row_cell_info_diamond['height'], self.col_overlay_var, self.row_overlay_var)
		self.cuts_man.in_central_rect_region[percentage] = '(' + self.gridAreas.center_rectangles[percentage]['name'] + ')'
		self.cuts_man.out_central_rect_region[percentage] = '(!' + self.gridAreas.center_rectangles[percentage]['name'] + ')'

	def CreateGridText(self):
		"""
		Creates the text that can be placed in profile diamond plots to easily navigate rows and columns
		:return:
		"""
		lims_dic = self.GetDiamondMapPlotLimits()
		xmin, xmax, ymin, ymax = lims_dic['xmin'], lims_dic['xmax'], lims_dic['ymin'], lims_dic['ymax']
		binsx, binsy = int(lims_dic['binsx']), int(lims_dic['binsy'])
		miny = min(self.row_cell_info_diamond['0_even'], self.row_cell_info_diamond['0_odd'])
		self.gridTextDiamond = ro.TH2F('gridText_diamond', 'gridText_diamond', binsx, xmin, xmax, binsy, ymin, ymax)
		for col in xrange(0, self.num_cols):
			x0 = self.dia_cols.cols[col].cells[0].xcenter
			y0 = miny
			self.gridTextDiamond.Fill(x0, y0 - 0.1 * self.row_cell_info_diamond['height'], (col + 0.01))
		max_num_rows = max(self.row_cell_info_diamond['num_even'], self.row_cell_info_diamond['num_odd'])
		for row in xrange(0, max_num_rows):
			if row < len(self.dia_cols.cols[0].cells):
				x0 = self.dia_cols.cols[0].cells[row].xcenter - self.dia_cols.cols[0].cells[row].w / 2.0
				y0 = self.dia_cols.cols[0].cells[row].ycenter
				self.gridTextDiamond.Fill(x0 - 0.1, y0, (row + 0.01))
			if self.num_cols > 1:
				if row < len(self.dia_cols.cols[1].cells):
					y0 = self.dia_cols.cols[1].cells[row].ycenter
					self.gridTextDiamond.Fill(x0 - 0.6, y0, (row + 0.01))
		self.gridTextDiamond.SetMarkerSize(0.8)

	def GetDiamondMapPlotLimits(self, border=3):
		"""
		Calculates the limits for profile 2D Diamond plots
		:param border: number of cells to leave as a border for the plot
		:return: dictionary with the calculated parameters
		"""
		xmin = self.dia_cols.cols[0].xcenter - self.row_cell_info_diamond['width'] * (border + 0.5) / self.col_pitch
		xmax = self.dia_cols.cols[-1].xcenter + self.row_cell_info_diamond['width'] * (border + 0.5) / self.col_pitch
		binsx = RoundInt((xmax - xmin) * self.bins_per_ch_x)
		ymin = min(self.row_cell_info_diamond['0_even'], self.row_cell_info_diamond['0_odd']) + self.row_cell_info_diamond['y_off'] - self.row_cell_info_diamond['height'] * border
		ymax = max(self.row_cell_info_diamond['up_even'], self.row_cell_info_diamond['up_odd']) + self.row_cell_info_diamond['y_off'] + self.row_cell_info_diamond['height'] * border
		binsy = RoundInt((ymax - ymin) * self.bins_per_ch_y / self.row_cell_info_diamond['height'])
		return {'xmin': xmin, 'xmax': xmax, 'binsx': binsx, 'ymin': ymin, 'ymax': ymax, 'binsy': binsy}

	def DrawProfile2DDiamond(self, name, varz='clusterChargeN', varname='PH [ADC]', cells='', cuts='', transp_ev=True, plot_option='prof colz'):
		"""
		Method used to make a profile map fo the diamond using the self.bins_per_ch_{x/y} variables to indicate how many bins are used for each column or row
		:param name: name of the profile. This grants easy access of the profile from the dictionary self.profile
		:param varz: z-variable used for the profile 2D which normally is a PH variable
		:param varname: Name to be displayed in the Color bar of the plot
		:param cells: indicates which cut should be applied for the fiducial region. 'good' is the selected cells, 'bad' is the non selected cells, 'all' are all the cells in the transparent grid,
				and '' is show all the region including the borders of the diamond
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		lims_dic = self.GetDiamondMapPlotLimits()
		xmin, xmax, ymin, ymax = lims_dic['xmin'], lims_dic['xmax'], lims_dic['ymin'], lims_dic['ymax']
		deltax, deltay = float(xmax - xmin) / lims_dic['binsx'], float(ymax - ymin) / lims_dic['binsy']
		xname = 'dia X ch'
		yname = 'sil pred Y [#mum]'
		tempc = self.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
		self.DrawProfile2D(name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, 'diaChXPred', 'diaChYPred', varz, varname, tempc, transp_ev, plot_option)

	def DrawProfile2DDiamondChannel(self, name, varx='clusterChannel0', xname='C0', varz='clusterChargeN', varname='PH [ADC]', cells='', cuts='', transp_ev=True, plot_option='prof colz'):
		"""
		Method used to plot the profile of a whole channel, to help identify defective cells, or clustered negative contributions
		:param name: name of the profile plot
		:param varx: variable to use in the X axis
		:param xname: name of the variable in the x axis
		:param varz: variable in the color palette. Most of the times it is a PH variable
		:param varname: name of the variable in z
		:param cells: indicates which cut should be applied for the fiducial region. 'good' is the selected cells, 'bad' is the non selected cells, 'all' are all the cells in the transparent grid,
				and '' is show all the region including the borders of the diamond
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param transp_ev:  Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		lims_dic = self.GetDiamondMapPlotLimits()
		xmin, xmax, ymin, ymax = lims_dic['xmin'], lims_dic['xmax'], lims_dic['ymin'], lims_dic['ymax']
		deltax, deltay = 1, float(ymax - ymin) / lims_dic['binsy']
		yname = 'sil pred Y [#mum]'
		tempc = self.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
		self.DrawProfile2D(name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, 'diaChYPred', varz, varname, tempc, transp_ev, plot_option)

	def DrawProfile1D(self, name, xmin, xmax, deltax, xname, varx, vary='clusterChargeN', yname='PH[ADC]', cuts='', transp_ev=True, plot_option='prof e hist'):
		"""
		Method to draw 1D Profiles
		:param name: name of the profile
		:param xmin: min x
		:param xmax: max x
		:param deltax: delta x
		:param xname: name of the x axis
		:param varx: variable used in the x axis
		:param vary: variable used in the y axis
		:param yname: name of the y axis
		:param cuts: string with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof e hist creates a profile plot with error bars
		:return:
		"""
		ro.TFormula.SetMaxima(100000)
		if name in self.profile.keys():
			self.profile[name].Delete()

		self.profile[name] = ro.TProfile('h_' + name, 'h_' + name, int(RoundInt((xmax - xmin)/float(deltax), 'f8') + 2), xmin - deltax, xmax + deltax)
		self.profile[name].GetXaxis().SetTitle(xname)
		self.profile[name].GetYaxis().SetTitle(yname)
		if 'goff' not in plot_option:
			if name in self.canvas.keys():
				self.canvas[name].Close()
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		list_cuts = ['transparentEvent'] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.trans_tree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), temp_cut, plot_option)
		ro.gPad.Update()
		SetDefault1DStats(self.profile[name])
		ro.TFormula.SetMaxima(1000)

	def DrawProfile2D(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, varz='clusterChargeN', zname='PH[ADC]', cuts='', transp_ev=True, plot_option='colz prof'):
		"""
		Method used by many other plotting methods to create 2D profiles
		:param name: name of the profile plot
		:param xmin: min x
		:param xmax: max x
		:param deltax: delta x
		:param xname: name of x axis
		:param ymin: min y
		:param ymax: max y
		:param deltay: delta y
		:param yname: name of y axis
		:param varx: variable used in the x axis
		:param vary: variable used in the y axis
		:param varz: variable used in the z axis
		:param zname: name of the z axis
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		ro.TFormula.SetMaxima(100000)
		if name in self.profile.keys():
			self.profile[name].Delete()

		self.profile[name] = ro.TProfile2D('h_' + name, 'h_' + name, int(RoundInt((xmax - xmin)/float(deltax), 'f8') + 2), xmin - deltax, xmax + deltax, int(RoundInt((ymax - ymin)/float(deltay), 'f8') + 2), ymin - deltay, ymax + deltay)
		self.profile[name].GetXaxis().SetTitle(xname)
		self.profile[name].GetYaxis().SetTitle(yname)
		self.profile[name].GetZaxis().SetTitle(zname)
		if 'goff' not in plot_option:
			if name in self.canvas.keys():
				self.canvas[name].Close()
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		list_cuts = [self.cuts_man.transp_ev] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.trans_tree.Draw('{z}:{y}:{x}>>h_{n}'.format(z=varz, y=vary, x=varx, n=name), temp_cut, plot_option)
		if 'goff' not in plot_option:
			ro.gPad.Update()
		SetDefault2DStats(self.profile[name])
		ro.TFormula.SetMaxima(1000)

	def DrawTCutGs(self, name, typ):
		"""
		Method to plot the TCutGs that delimit the cells in the transparent grid.
		:param name: name of the histogram or profile where to plot the grid
		:param typ: if typ is 'diamond' it will plot all the cutg objects on each cell and the grid text. If typ is 'centers' it will only plot the cutg_center of all the cells
		:return:
		"""
		self.canvas[name].cd()
		ro.gStyle.SetPaintTextFormat(".0f")
		if typ == 'diamond':
			self.gridTextDiamond.Draw('same TEXT0')
		if name in self.profile.keys():
			self.profile[name].Draw('same colz')
		elif name in self.histo.keys():
			self.histo[name].Draw('same colz')
		for col in xrange(self.num_cols):
			for row in xrange(self.dia_cols.cols[col].num_rows):
				if typ == 'diamond':
					self.dia_cols.cols[col].cells[row].cutg.Draw('same')
				elif typ == 'centers':
					self.dia_cols.cols[col].cells[row].cutg_center.Draw('same')
		ro.gPad.Update()

	def GetOccupancyFromProfile(self, name, plot_option='colz'):
		"""
		Creates a 2D histogram from a given profile of the occupancy of the 2D profile
		:param name: name of the targe profile 2D plot
		:param plot_option: plotting option used for root. colz creates a color bar for a Histo 2D plot
		:return:
		"""
		name_occupancy = 'hit_map_' + name
		if name_occupancy in self.histo.keys():
			self.histo[name_occupancy].Delete()
		self.histo[name_occupancy] = self.profile[name].ProjectionXY('h_' + name_occupancy, 'B')
		self.histo[name_occupancy].SetTitle('h_' + name_occupancy)
		self.histo[name_occupancy].GetXaxis().SetTitle(self.profile[name].GetXaxis().GetTitle())
		self.histo[name_occupancy].GetYaxis().SetTitle(self.profile[name].GetYaxis().GetTitle())
		self.histo[name_occupancy].GetZaxis().SetTitle('entries')
		if 'goff' not in plot_option:
			if name_occupancy in self.canvas.keys():
				self.canvas[name_occupancy].Close()
			self.canvas[name_occupancy] = ro.TCanvas('c_' + name_occupancy, 'c_' + name_occupancy, 1)
			self.canvas[name_occupancy].cd()
			self.histo[name_occupancy].Draw(plot_option)
			ro.gPad.Update()
			SetDefault2DStats(self.histo[name_occupancy])

	def DrawHisto2D(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, cuts='', transp_ev=True):
		"""
		Method used by other plotting methods to generate 2D histograms
		:param name: name of the histogram
		:param xmin: min x
		:param xmax: max x
		:param deltax: delta x
		:param xname: name of x axis
		:param ymin: min y
		:param ymax: max y
		:param deltay: delta y
		:param yname: name of y axis
		:param varx: variable used for the x axis
		:param vary: variable used for the y axis
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:return:
		"""
		ro.TFormula.SetMaxima(100000)
		if name in self.histo.keys():
			self.histo[name].Delete()
		self.histo[name] = ro.TH2F('h_' + name, 'h_' + name, int(RoundInt((xmax - xmin) / deltax, 'f8') + 2), xmin - deltax, xmax + deltax, int(RoundInt((ymax - ymin) / deltay, 'f8') + 2), ymin - deltay, ymax + deltay)
		self.histo[name].GetXaxis().SetTitle(xname)
		self.histo[name].GetYaxis().SetTitle(yname)
		self.histo[name].GetZaxis().SetTitle('entries')
		if name in self.canvas.keys():
			self.canvas[name].Close()
		self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		self.canvas[name].cd()
		list_cuts = [self.cuts_man.transp_ev] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cut = '&&'.join(list_cuts)
		self.trans_tree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), temp_cut, 'colz')
		ro.gPad.Update()
		SetDefault2DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def DrawHisto1D(self, name, xmin, xmax, deltax, var='clusterChargeN', varname='PH[ADC]', cuts='', transp_ev=True, option='e hist'):
		"""
		Method used by other plotting methods to generate 1D histograms
		:param name: name of the histogram
		:param xmin: min x
		:param xmax: max x
		:param deltax: delta x
		:param var: variable used for the histogram
		:param varname: name of the x axis
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param option: plotting option used for root. e hist creates a histogram bar plot with error bars
		:return:
		"""
		ro.TFormula.SetMaxima(100000)
		if name in self.histo.keys():
			self.histo[name].Delete()
		self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, int(RoundInt((xmax - xmin) / float(deltax))), xmin, xmax)
		self.histo[name].GetXaxis().SetTitle(varname)
		self.histo[name].GetYaxis().SetTitle('entries')
		list_cuts = ['transparentEvent'] if transp_ev else []
		if cuts != '':
			list_cuts.append(cuts)
		temp_cuts = '&&'.join(list_cuts)
		if 'goff' not in option:
			if name in self.canvas.keys():
				self.canvas[name].Close()
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		self.trans_tree.Draw('{v}>>h_{n}'.format(v=var, n=name), temp_cuts, option)
		if 'goff' not in option:
			SetDefault1DCanvasSettings(self.canvas[name])
			ro.gPad.Update()
			SetDefault1DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def GetMeanPHPerCell(self, var='clusterChargeN', typ='adc'):
		"""
		This method tries to load a pickle with the average value of var of type typ ('adc' or 'snr'). If the pickle does not have that info, it will be calculated and saved.
		The result is stored in the dictionary self.mean_ph_cell_dic[typ]
		:param var: variable to calculate the mean ph of each cell
		:param typ: it tells if it is in snr or adc. It can take values 'adc' or 'snr'
		:return:
		"""
		means_temp = {}
		if os.path.isfile('{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir, t=typ)):
			with open('{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir, t=typ)) as pkl:
				means_temp = pickle.load(pkl)
				if means_temp:
					if var in means_temp.keys():
						self.mean_ph_cell_dic[typ] = means_temp[var]
						print 'Loaded pickle with mean PH info for each cell'
						return
					else:
						print 'The existing file does not have the info for {v2}. Will calculate!'.format(v2=var)
					# if 'var' in means_temp.keys():
					# 	if means_temp['var'] == var:
					# 		self.mean_ph_cell_dic[typ] = means_temp['dic']
					# 		print 'Loaded pickle with mean PH info for each cell'
					# 		return
					# 	else:
					# 		print 'The existing file has info for {v} but the requested variable is {v2}. Will recalculate!'.format(v=means_temp['var'], v2=var)
		print 'Calculating the mean PH value for each cell:'
		numcells = int(self.dia_cols.num_cols * max(self.row_cell_info_diamond['num_even'], self.row_cell_info_diamond['num_odd']))
		tempbar = CreateProgressBarUtils(numcells)
		tempbar.start()
		for col in xrange(self.dia_cols.num_cols):
			self.mean_ph_cell_dic[typ][col] = {}
			for row in xrange(self.dia_cols.cols[col].num_rows):
				if typ == 'adc':
					self.trans_tree.Draw(var + '>>temphrc(200,0,4800)', 'transparentEvent&&({n})'.format(n=self.dia_cols.cols[col].cells[row].cutg.GetName()), 'goff')
				else:
					self.trans_tree.Draw(var + '>>temphrc(200,0,480)', 'transparentEvent&&({n})'.format(n=self.dia_cols.cols[col].cells[row].cutg.GetName()), 'goff')
				temph = ro.gDirectory.Get('temphrc')
				self.mean_ph_cell_dic[typ][col][row] = temph.GetMean()
				temph.Reset('ICES')
				temph.Delete()
				del temph
				tempbar.update(col * self.dia_cols.cols[col].num_rows + row + 1)
		tempbar.finish()
		meanph_obj = means_temp
		meanph_obj[var] = self.mean_ph_cell_dic[typ]
		if not os.path.isdir('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
		print 'Saving limits in pickle file', '{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir, t=typ), '...', ; sys.stdout.flush()
		pickle.dump(meanph_obj, open('{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir, t=typ), 'wb'))
		print 'Done'

	def SelectGoodAndBadByThreshold(self, val=500, var='clusterChargeN', typ='adc'):
		"""
		This method will classify each cell as 'good' or 'bad' (selection or not in selection) depending on a threshold given by val, a variable var of type typ.
		The results are stored in the object gridAreas
		:param val: threshold to use for classification
		:param var: variable of ph used for discrimination
		:param typ: type 'adc' or 'snr' used for classification. The variable var has to be of this type.
		:return:
		"""
		if len(self.mean_ph_cell_dic[typ].keys()) > 0:
			for col in xrange(self.dia_cols.num_cols):
				for row in xrange(self.dia_cols.cols[col].num_rows):
					if self.mean_ph_cell_dic[typ][col][row] > val:
						self.gridAreas.AddGoodAreas(col, row, self.dia_cols)
					else:
						self.gridAreas.AddBadAreas(col, row, self.dia_cols)
		else:
			self.GetMeanPHPerCell(var, typ)
			self.SelectGoodAndBadByThreshold(val, var)

	def FindThresholdCutFromCells(self, var='clusterChargeN', typ='adc', xmin=0, xmax=4800, deltax=10, it=0):
		"""
		This method finds the threshold cut using the distribution of PH for each cell. The threshold is set as a number of sigmas (self.threshold_criteria_in_sigmas) below the mean of the gaussian fit of the distribution of PH
		:param var: variable of type 'typ' used for classification
		:param typ: type of variable 'var' used. It has to be either 'adc' or 'snr'
		:param xmin: min value of the distribution of PH determined by var
		:param xmax: max value of the distribution of PH determined by var
		:param deltax: bin width of the distribution of PH determined by var
		:param it: number of iterations in this method to prevent infinite recursion
		:return:
		"""
		if len(self.mean_ph_cell_dic[typ].keys()) > 0:
			self.DrawMeanPHCellsHisto(var, typ, xmin, xmax, deltax)
			self.FitGaus('mean_ph_per_cell_' + typ, 2, 5)
			if GetHistoAverageBinContent(self.histo['mean_ph_per_cell_' + typ]) < 10 and it < 20:
				self.FindThresholdCutFromCells(var, typ, xmin, xmax, deltax * 1.25, it + 1)
			elif self.fits['mean_ph_per_cell_' + typ].Ndf() < 2 and it < 20:
				self.FindThresholdCutFromCells(var, typ, xmin, xmax, deltax * 0.8, it + 1)
			elif 0.5 <= self.fits['mean_ph_per_cell_' + typ].Chi2() / self.fits['mean_ph_per_cell_' + typ].Ndf() < 0.9 and it < 20:
				self.FindThresholdCutFromCells(var, typ, xmin, xmax, deltax * 0.8, it + 1)
		else:
			self.GetMeanPHPerCell(var, typ)
			self.FindThresholdCutFromCells(var, typ, xmin, xmax, deltax, it + 1)
		if typ == 'adc':
			self.threshold = self.fits['mean_ph_per_cell_' + typ].Parameter(1) - self.threshold_criteria_in_sigmas * self.fits['mean_ph_per_cell_' + typ].Parameter(2)
			if self.threshold <= deltax:
				self.threshold_criteria_in_sigmas /= 2.0
				self.FindThresholdCutFromCells(var, typ, xmin, xmax, deltax, it)
		else:
			self.threshold_snr = self.fits['mean_ph_per_cell_' + typ].Parameter(1) - self.threshold_criteria_in_sigmas * self.fits['mean_ph_per_cell_' + typ].Parameter(2)
		thresh = self.threshold if typ == 'adc' else self.threshold_snr
		self.line['threshold_' + typ] = ro.TLine(thresh, 0, thresh, self.histo['mean_ph_per_cell_' + typ].GetMaximum())
		self.line['threshold_' + typ].SetLineColor(ro.kBlue)
		self.line['threshold_' + typ].SetLineWidth(2)
		self.line['threshold_' + typ].SetLineStyle(2)
		self.line['threshold_' + typ].Draw('same')

	def DrawMeanPHCellsHisto(self, var='clusterChargeN', typ='adc', xmin=0, xmax=4800, deltax=50, draw_opt='e hist', supress0=True):
		"""
		This method plots the PH distribution of variable 'var' of type 'typ' with at least 10 bins
		:param var: variable of type 'typ' used for the distribution
		:param typ: type of variable 'var'. Has to be either 'adc' or 'snr'
		:param xmin: min value of the distribution of PH determined by var
		:param xmax: max value of the distribution of PH determined by var
		:param deltax: bin width of the distribution of PH determined by var
		:param draw_opt: draw option of the histogram. Default 'e hist'
		:return:
		"""
		if len(self.mean_ph_cell_dic[typ].keys()) > 0:
			nameh = 'mean_ph_per_cell_' + typ
			if nameh in self.histo.keys():
				self.histo[nameh].Reset('ICES')
				self.histo[nameh].Delete()
				del self.histo[nameh]
			limsh = Get1DLimits(xmin, xmax, deltax) if not supress0 else Get1DLimits(deltax, xmax, deltax)
			self.histo[nameh] = ro.TH1F('h_' + nameh, 'h_' + nameh, int(RoundInt((limsh['max'] - limsh['min']) / deltax)), limsh['min'], limsh['max'])
			for col, rowMean in self.mean_ph_cell_dic[typ].iteritems():
				for row, meanh in rowMean.iteritems():
					self.histo[nameh].Fill(meanh)
			if 'goff' not in draw_opt.lower():
				if nameh not in self.canvas.keys():
					self.canvas[nameh] = ro.TCanvas('c_' + nameh, 'c_' + nameh, 1)
				self.canvas[nameh].cd()
				self.histo[nameh].Draw(draw_opt)
				SetDefault1DCanvasSettings(self.canvas[nameh])
				ro.gPad.Update()
			SetDefault1DStats(self.histo[nameh])
			if(self.histo[nameh].FindLastBinAbove() - self.histo[nameh].FindFirstBinAbove() + 1) < 10:
				self.DrawMeanPHCellsHisto(var, typ, xmin, xmax, deltax * 0.8, draw_opt, supress0)
		else:
			self.GetMeanPHPerCell(var, typ)
			self.DrawMeanPHCellsHisto(var, typ, xmin, xmax, deltax, draw_opt, supress0)

	def DrawCentralArea(self, name, percent):
		"""
		Method to draw the central rectangle in the percent studies
		:param name: name of the canvas where the central region will be drawn
		:param percent: percentage of the inner area inside the rectangle to be plotted
		:return:
		"""
		self.canvas[name].cd()
		if percent in self.gridAreas.center_rectangles.keys():
			self.gridAreas.center_rectangles[percent]['tcutg'].Draw('same')

	def DrawGoodAreasDiamondCenters(self, name):
		"""
		Method used to plot the central cutg of each cell in the selected regions ('good'). It is used in general to highlight the selected cells
		:param name: name of the plot where the central regions will be drawn
		:return:
		"""
		self.DrawGoodAreas(name, typ='centers')

	def DrawGoodAreas(self, name, typ):
		"""
		Method used to plot cutgs of goodAreas or goodAreas centers on a specified plot
		:param name: name of the plot where the cutgs will be drawn
		:param typ: if typ is 'centers' then it will only plot the cutg_center of the good cells (selected). if typ is 'diamond', it will plat the cutg that delimits each good cell (selected region)
		:return:
		"""
		self.canvas[name].cd()
		if typ == 'diamond':
			for area in self.gridAreas.goodAreas_diamond:
				area.Draw('same')
		elif typ == 'centers':
			for area in self.gridAreas.goodAreas_diamond_centers:
				area.Draw('same')

	def DrawBadAreasDiamond(self, name):
		"""
		Method used to plot the cutg of each cell in the non-selected regions ('bad')
		:param name: name of the plot where the regions will be drawn
		:return:
		"""
		self.DrawBadAreas(name, type='diamond')

	def DrawBadAreasDiamondCenters(self, name):
		"""
		Method used to plot the central cutg of each cell in the non-selected regions ('bad'). It is used in general to highlight the non-selected cells
		:param name: name of the plot where the central regions will be drawn
		:return:
		"""
		self.DrawBadAreas(name, type='centers')

	def DrawBadAreas(self, name, type):
		"""
		Method used to plot cutgs of badAreas or badAreas centers on a specified plot
		:param name: name of the plot where the cutgs will be drawn
		:param type: if typ is 'centers' then it will only plot the cutg_center of the bad cells (non-selected). if typ is 'diamond', it will plat the cutg that delimits each bad cell (non selected region)
		:return:
		"""
		self.canvas[name].cd()
		if type == 'diamond':
			for area in self.gridAreas.badAreas_diamond:
				area.Draw('same')
		elif type == 'centers':
			for area in self.gridAreas.badAreas_diamond_centers:
				area.Draw('same')

	def ResetHistos(self):
		"""
		Clearse all the histograms stored in self.histo
		:return:
		"""
		for histo in self.histo.itervalues():
			histo.Delete()
		self.histo = {}

	def ResetProfiles(self):
		"""
		Clears all the profiles stored in self.profile
		:return:
		"""
		for profile in self.profile.itervalues():
			profile.Delete()
		self.profile = {}

	def ResetCanvas(self):
		"""
		Clears all the canvas stored in self.canvas
		:return:
		"""
		for canvas in self.canvas.itervalues():
			canvas.Clear()
			canvas.Close()
		self.canvas = {}

	def ResetPlots(self):
		"""
		Clears all the plots done!
		:return:
		"""
		self.ResetHistos()
		self.ResetProfiles()
		self.ResetCanvas()

	def AddGoodAreas(self, col, row):
		"""
		Tags a cell identified by column 'col' and row 'row' as good (in selection)
		:param col: cell's column number
		:param row: cell's row number
		:return:
		"""
		self.gridAreas.AddGoodAreas(col, row, self.dia_cols)

	def AddBadAreas(self, col, row):
		"""
		Tags a cell identified by column 'col' and row 'row' as bad (not in selection)
		:param col: cell's column number
		:param row: cell's row number
		:return:
		"""
		self.gridAreas.AddBadAreas(col, row, self.dia_cols)

	def AddGoodAreasRow(self, row, coli=0, colf=0):
		"""
		Tags all the cells between 'coli' and 'colf' inclusively in the row 'row' as good (in selection)
		:param row: number of the row used for the addition
		:param coli: starting column
		:param colf: finishing column
		:return:
		"""
		self.gridAreas.AddGoodAreasRow(row, coli, colf, self.dia_cols)

	def AddGoodAreasCol(self, col, rowi=0, rowf=0):
		"""
		Tags all the cells between 'rowi' and 'rowf' inclusively in the column 'col' as good (in selection)
		:param col: number of the column used for the addition
		:param rowi: starting row
		:param rowf: finishing row
		:return:
		"""
		self.gridAreas.AddGoodAreasCol(col, rowi, rowf, self.dia_cols)

	def AddRemainingToBadAreas(self):
		"""
		Calls the gridAreas method to mark the remaining cells as bad (not in selection)
		:return:
		"""
		self.gridAreas.AddRemainingToBadAreas(self.dia_cols)

	def RemoveFromGoodArea(self, col, row):
		"""
		Tags a cell identified by column 'col' and row 'row' as bad (not in selection) if it was previously tagged as good (in selection) it is removed from the corresponding lists
		:param col: cell's column nomber
		:param row: cell's row number
		:return:
		"""
		self.gridAreas.RemoveFromGoodArea(col, row, self.dia_cols)

	def ResetAreas(self):
		"""
		Calls the method from gridAreas object to reset the areas
		:return:
		"""
		self.gridAreas.ResetAreas()

	def DrawPHInArea(self, name, var='clusterChargeN', cells='all', cuts='', varname='PH[ADC]', xmin=10000, xmax=-10000, deltax=-1, typ='adc', drawHisto=True):
		"""
		Method used Draw the PH histogram in a certain region aka area (good, bad, all)
		:param name: name of the histogram
		:param var: variable used to create the histogram. In general it is a PH variable
		:param cells: the options are 'all', 'good', 'bad' to select the region where to plot the histogram
		:param cuts: srting with the cuts to be applied by root
		:param varname: name of the X axis
		:param xmin: min x. if xmin is 10000, it will use the default phmin
		:param xmax: max x. if xmax is -10000, it will use the default phmax
		:param deltax: delta x. If deltax is -1, it will used the default phbins to calculate the deltax
		:param typ: if indicates if 'adc' is requested, or if 'snr' is requested
		:param drawHisto: if true, the histogram will be drawn in a canvas
		:return:
		"""
		if cells == 'good':
			self.DrawPHGoodAreas(name, var, cuts, varname=varname, xmin=xmin, xmax=xmax, deltax=deltax, typ=typ, drawHisto=drawHisto)
		elif cells == 'bad':
			self.DrawPHBadAreas(name, var, cuts, varname=varname, xmin=xmin, xmax=xmax, deltax=deltax, typ=typ, drawHisto=drawHisto)
		else:
			self.DrawPHAllAreas(name, var, cuts, varname=varname, xmin=xmin, xmax=xmax, deltax=deltax, typ=typ, drawHisto=drawHisto)

	def DrawPHGoodAreas(self, name, var='clusterChargeN', cuts='', transp_ev=True, varname='PH[ADC]', xmin=10000, xmax=-10000, deltax=-1, typ='adc', drawHisto=True):
		"""
		Method to draw the PH histogram in the cells marked as good
		:param name: name of the histogram
		:param var: variable used to create the histogram. In general it is a PH variable
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param varname: name of the X axis
		:param xmin: min x. if xmin is 10000, it will use the default phmin
		:param xmax: max x. if xmax is -10000, it will use the default phmax
		:param deltax: delta x. If deltax is -1, it will used the default phbins to calculate the deltax
		:param typ: if indicates if 'adc' is requested, or if 'snr' is requested
		:param drawHisto: if true, the histogram will be drawn in a canvas
		:return:
		"""
		temp_cut = self.cuts_man.ConcatenateCutWithCells(cuts, 'good')
		phmin = xmin if xmin != 10000 else self.phmin if typ == 'adc' else self.phmin / 10.
		phmax = xmax if xmax != -10000 else self.phmax if typ == 'adc' else self.phmax / 10.
		deltx = deltax if deltax != -1 else float(phmax - phmin) / float(self.phbins)
		graphopt = 'e hist' if drawHisto else 'e hist goff'
		self.DrawHisto1D(name, phmin, phmax, deltx, var, varname, temp_cut, transp_ev, option=graphopt)

	def DrawPHAllAreas(self, name, var='clusterChargeN', cuts='', transp_ev=True, varname='PH[ADC]', xmin=10000, xmax=-10000, deltax=-1, typ='adc', drawHisto=True):
		"""
		Method to draw the PH histogram in all the cells in the transparent grid
		:param name: name of the histogram
		:param var: variable used to create the histogram. In general it is a PH variable
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param varname: name of the X axis
		:param xmin: min x. if xmin is 10000, it will use the default phmin
		:param xmax: max x. if xmax is -10000, it will use the default phmax
		:param deltax: delta x. If deltax is -1, it will used the default phbins to calculate the deltax
		:param typ: if indicates if 'adc' is requested, or if 'snr' is requested
		:param drawHisto: if true, the histogram will be drawn in a canvas
		:return:
		"""
		temp_cut = self.cuts_man.ConcatenateCutWithCells(cuts, 'all')
		phmin = xmin if xmin != 10000 else self.phmin if typ == 'adc' else self.phmin / 10.
		phmax = xmax if xmax != -10000 else self.phmax if typ == 'adc' else self.phmax / 10.
		deltx = deltax if deltax != -1 else float(phmax - phmin) / float(self.phbins)
		graphopt = 'e hist' if drawHisto else 'e hist goff'
		self.DrawHisto1D(name, phmin, phmax, deltx, var, varname, temp_cut, transp_ev, option=graphopt)

	def DrawPHBadAreas(self, name, var='clusterChargeN', cuts='', transp_ev=True, varname='PH[ADC]', xmin=10000, xmax=-10000, deltax=-1, typ='adc', drawHisto=True):
		"""
		Method to draw the PH histogram in all the cells marked as bad
		:param name: name of the histogram
		:param var: variable used to create the histogram. In general it is a PH variable
		:param cuts: srting with the cuts to be applied by root
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param varname: name of the X axis
		:param xmin: min x. if xmin is 10000, it will use the default phmin
		:param xmax: max x. if xmax is -10000, it will use the default phmax
		:param deltax: delta x. If deltax is -1, it will used the default phbins to calculate the deltax
		:param typ: if indicates if 'adc' is requested, or if 'snr' is requested
		:param drawHisto: if true, the histogram will be drawn in a canvas
		:return:
		"""
		temp_cut = self.cuts_man.ConcatenateCutWithCells(cuts, 'bad')
		phmin = xmin if xmin != 10000 else self.phmin if typ == 'adc' else self.phmin / 10.
		phmax = xmax if xmax != -10000 else self.phmax if typ == 'adc' else self.phmax / 10.
		deltx = deltax if deltax != -1 else float(phmax - phmin) / float(self.phbins)
		graphopt = 'e hist' if drawHisto else 'e hist goff'
		self.DrawHisto1D(name, phmin, phmax, deltx, var, varname, temp_cut, transp_ev, option=graphopt)

	def DrawProfile2DDiamondChannelOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True, plot_option='prof colz'):
		"""
		Method used to plot 2D profile of overlaid columns (channels).
		:param name: name of the overlaid 2D profile
		:param var: z-variable used for the profile 2D which normally is a PH variable
		:param cells: indicates which cut should be applied for the fiducial region. 'good' is the selected cells, 'bad' is the non selected cells, 'all' are all the cells in the transparent grid,
				and '' is show all the region including the borders of the diamond
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		temp_cuts = self.cuts_man.ConcatenateCutWithCells(cuts, cells)
		lims_dic = self.GetDiamondMapPlotLimits()
		ymin, ymax = lims_dic['ymin'], lims_dic['ymax']
		deltay = float(ymax - ymin) / lims_dic['binsy']
		xmin, xmax, deltax = -self.row_cell_info_diamond['width'] / 2.0, self.row_cell_info_diamond['width'] / 2.0, self.cell_resolution
		self.DrawProfile2D(name, xmin, xmax, deltax, 'dia X [#mum]', ymin, ymax, deltay, 'dia Y [#mum]', self.col_overlay_var2, 'diaChYPred', var, 'PH [ADC]', temp_cuts, transp_ev, plot_option)

	def DrawProfile2DDiamondRowOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True, plot_option='prof colz'):
		"""
		Method used to plot 2D profile of overlaid rows.
		:param name: name of the overlaid 2D profile
		:param var: z-variable used for the profile 2D which normally is a PH variable
		:param cells: indicates which cut should be applied for the fiducial region. 'good' is the selected cells, 'bad' is the non selected cells, 'all' are all the cells in the transparent grid,
				and '' is show all the region including the borders of the diamond
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		lims_dic = self.GetDiamondMapPlotLimits()
		xmin, xmax = lims_dic['xmin'], lims_dic['xmax']
		deltax = float(xmax - xmin) / lims_dic['binsx']
		ymin, ymax, deltay = -self.row_cell_info_diamond['height'] / 2.0, self.row_cell_info_diamond['height'] / 2.0, self.cell_resolution
		temp_cuts = self.cuts_man.ConcatenateCutWithCells(cuts, cells)
		self.DrawProfile2D(name, xmin, xmax, deltax, 'dia X ch', ymin, ymax, deltay, 'dia Y [#mum]', 'diaChXPred', self.row_overlay_var2, var, 'PH [ADC]', temp_cuts, transp_ev, plot_option)

	def DrawProfile2DDiamondCellOverlay(self, name, var='clusterChargeN', cells='all', cuts='', transp_ev=True, plot_option='prof colz', varname='PH [ADC]', typ='adc'):
		"""
		Method used to plot 2D profile of overlaid cells.
		:param name: name of the overlaid 2D profile
		:param var: z-variable used for the profile 2D which normally is a PH variable
		:param cells: indicates which cut should be applied for the fiducial region. 'good' is the selected cells, 'bad' is the non selected cells, 'all' are all the cells in the transparent grid,
				and '' is show all the region including the borders of the diamond
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param transp_ev: Most of the times it is true. Sets the trnasparentEvent cut. Events which don't have the transparentEvent cut, normally are useless
		:param plot_option: plotting option used for root. prof colz creates a color bar for a profile 2D plot
		:return:
		"""
		xmin, xmax, deltax = -self.row_cell_info_diamond['width'] / 2.0, self.row_cell_info_diamond['width'] / 2.0, self.cell_resolution
		ymin, ymax, deltay = -self.row_cell_info_diamond['height'] / 2.0, self.row_cell_info_diamond['height'] / 2.0, self.cell_resolution
		temp_cuts = self.cuts_man.ConcatenateCutWithCells(cuts, cells)
		self.DrawProfile2D(name, xmin, xmax, deltax, 'dia X [#mum]', ymin, ymax, deltay, 'dia Y [#mum]', self.col_overlay_var2, self.row_overlay_var2, var, varname, temp_cuts, transp_ev, plot_option)
		if typ == 'adc':
			self.profile[name].GetZaxis().SetRangeUser(self.phmin, self.phmax)
		else:
			self.profile[name].GetZaxis().SetRangeUser(self.phmin / 10., self.phmax / 10.)
		if 'goff' not in plot_option:
			ro.gPad.Update()

	def DrawEfficiencyGraph(self, name, var, cells, cuts, xmin, xmax, deltax, typ='adc', ymin_plot=0, sigma_errbar=ro.TMath.Erf(1/np.sqrt(2)), subf=4000.0):
		"""
		This method plots an efficiency graph for a certain PH variable of either 'adc' type or 'snr' type. The error bars are calculated numericaly using Bayes theorem to generate assymetric error bars
		:param name: name of efficiency graph
		:param var: variable used for efficiency discrimination
		:param cells: cells used: 'good', 'bad', 'all' or ''
		:param cuts: is a string which contains the cuts used besides the fiducial cut
		:param xmin: min x
		:param xmax: max x
		:param deltax: delta x
		:param typ: type of variable used for discriminating the efficiency: either 'adc' or 'snr'
		:param ymin_plot: minimum value of efficiency to show. This will reduce the number of calculated efficiencies
		:param sigma_errbar: error bar used for limit values
		:param subf: subfactor used to calculate the subdivisions between efficiency steps if efficiency_subdiv == 1
		:return:
		"""
		cut = self.cuts_man.transp_ev if cuts == '' else self.cuts_man.AndCuts([self.cuts_man.transp_ev, cuts])
		minimum = min(0, max(GetMinimumFromTree(self.trans_tree, var, cut), -9999))
		denominator = float(self.GetEventsForThCut(var, minimum, cells, cut))
		if denominator == 0:
			print 'The denominator to calculate the efficiencies is 0! Please check!'
			return
		xvalues = np.arange(xmin, xmax - (xmax - xmin) * ymin_plot ** 4, deltax, 'float64')
		print 'Getting efficiencies...', ; sys.stdout.flush()
		numerator = np.array(map(self.GetEventsForThCut, [var for i in xrange(xvalues.size)], xvalues, [cells for i in xrange(xvalues.size)], [cut for i in xrange(xvalues.size)]), 'float64')
		efficiency = np.divide(numerator, denominator, dtype='float64')
		print 'Done'
		# efficiencyDic = {th: numerator[th] / denominator for th in xvalues}
		lim_one_side = np.subtract(1, np.power(np.subtract(1, sigma_errbar, dtype='float64'), np.divide(1, denominator + 1, dtype='float64'), dtype='float64'), dtype='float64')
		(xinf, xsup, yinf, ysup) = (xmin - 10, xmax + 10, 0 - lim_one_side, 1 + lim_one_side) if ymin_plot == 0 else (xmin - 10, xmax + 10, ymin_plot - lim_one_side, 1 + lim_one_side)
		if ymin_plot != 0:
			for it, value in enumerate(xvalues):
				if efficiency[it] <= ymin_plot:
					xsup = value - deltax / 2.0
					yinf = ymin_plot - lim_one_side
					break
		self.efficiency_subdiv = int(max(RoundInt(subf / denominator), 1)) if self.efficiency_subdiv == 1 else self.efficiency_subdiv
		cont = 0
		while ((1.0 / (self.efficiency_subdiv * denominator)) < 1e-5 or (1.0 / (self.efficiency_subdiv * denominator)) > 1e-3 )and cont < 1000000:
			if (1.0 / (self.efficiency_subdiv * denominator)) < 1e-5:
				self.efficiency_subdiv = RoundInt(self.efficiency_subdiv * 0.95) if self.efficiency_subdiv > 1 else self.efficiency_subdiv * 0.5
				cont += 1
			elif (1.0 / (self.efficiency_subdiv * denominator)) > 1e-3:
				self.efficiency_subdiv = RoundInt(self.efficiency_subdiv * 1.5)
				cont += 1
		print '{c}. Calculate uncertainties for efficiency plot... the step used is {d}...'.format(c=cont, d=1.0 / (self.efficiency_subdiv * denominator)), ; sys.stdout.flush()
		ySigmas = map(self.FindAsymmetricUncertaintiesWithDiscrete, numerator, [denominator for i in xrange(numerator.size)], [sigma_errbar for i in xrange(numerator.size)])
		print 'Done'
		yLowerSigmas = np.array([ysigma['lower'] for ysigma in ySigmas], 'float64')
		yUpperSigmas = np.array([ysigma['upper'] for ysigma in ySigmas], 'float64')
		self.graph[name] = ro.TGraphAsymmErrors(len(xvalues), xvalues, efficiency, np.ones(xvalues.size, 'float64') * deltax / 2.0, np.ones(xvalues.size, 'float64') * deltax / 2.0, yLowerSigmas, yUpperSigmas)
		self.graph[name].SetNameTitle('g_' + name, 'g_' + name)
		self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
		self.graph[name].SetMarkerStyle(ro.TAttMarker.kFullDotMedium)
		self.graph[name].SetMarkerColor(ro.kRed)
		self.graph[name].Draw('AP')
		xtitle = 'Threshold [ADC]' if typ == 'adc' else 'Threshold [SNR]'
		self.graph[name].GetXaxis().SetTitle(xtitle)
		self.graph[name].GetYaxis().SetTitle('Efficiency')
		self.graph[name].GetXaxis().SetRangeUser(xinf, xsup)
		self.graph[name].GetYaxis().SetRangeUser(yinf, ysup)
		ro.gPad.Update()
		self.canvas[name].SetGridx()
		self.canvas[name].SetGridy()
		self.canvas[name].SetTicky()
		ro.gPad.Update()

	def GetEventsForThCut(self, var, threshold, cells, cut):
		"""
		Method used to get the number of entries that surpass a certain threshold for a certain variable after applying a cut on certain region of cells
		:param var: variable used for discrimination
		:param threshold: threshold used on variable var
		:param cells: cells used for calculation: 'good', 'all', 'bad', ''
		:param cut: is a string which contains the cuts used besides the fiducial cut
		:return: the number of entries that surpass the threshold
		"""
		temp_cut = self.cuts_man.GetThCut(var=var, th=threshold, cells=cells, cuts=cut, op='>=')
		return GetNumberEntriesFromTree(self.trans_tree, temp_cut)

	def BetaDistIntegral(self, k, n, xinf, xsup):
		"""
		Analytic approximation to calculate asymmetric uncertainties using Bayes theorem
		:param k:
		:param n:
		:param xinf:
		:param xsup:
		:return:
		"""
		return ro.TMath.BetaIncomplete(xsup, k + 1, n - k + 1) - ro.TMath.BetaIncomplete(xinf, k + 1, n - k + 1)

	def LagrangeFcn(self, npar, par):
		"""
		Lagrange function to minimize for analytic calculation of asymmetric uncertainties using Bayes theorem
		:param npar:
		:param par:
		:return:
		"""
		a = par[0]
		b = par[1]
		lambd = par[2]
		k = par[3]
		n = par[4]
		sigm = par[5]
		# if k == 0:
		# 	return float(b - lambd * (1 - np.power(1 - b, n + 1, dtype='float64') - sigm))
		# elif n == k:
		# 	return float(a - lambd * (1 - np.power(1 - a, n + 1, dtype='float64') - sigm))
		# else:
		return float(b + a - lambd * (self.BetaDistIntegral(k, n, float(k)/float(n) - a, float(k)/float(n) + b) - sigm))

	def MinuitFcn(self, npar, deriv, f, apar, iflag):
		""" meaning of parametrs:
		npar: number of parameters
		deriv: aray of derivatives df/dp_i (x), optional
		f: value of function to be minimised (typically chi2 or negLogL)
		apar: the array of parameters
		iflag: internal flag: 1 at first call , 3 at the last , 4 during minimisation
		"""
		f[0] = self.LagrangeFcn(npar, apar)

	def FindAsymmetricUncertaintiesWithMinuit(self, k, n, sigm, max_iter=100000, tolerance=0.1, minimizer='SIMPLEX', a0=0, b0=0, is_last=True):
		"""
		Attempt to use lagrange multipliers and minuit to estimate the asymmetric uncertainties. It was too slow and many times the solution did not converge. Discrete numerical approach is better
		:param k:
		:param n:
		:param sigm:
		:param max_iter:
		:param tolerance:
		:param minimizer:
		:param a0:
		:param b0:
		:param is_last:
		:return:
		"""
		myMinuit = ro.TMinuit(6)  # 6 parameters: a, b, lambd, k, n, sigma
		myMinuit.SetFCN(self.MinuitFcn)
		ro.gMinuit.Command('SET PRINT -1')
		myMinuit.SetPrintLevel(-1)
		ierflg = ro.Long(0)
		# Set Parameters
		initialValues = {}
		limSup = {}
		limInf = {}
		lim_binomial = np.sqrt(k * (1 - float(k) / float(n))) / float(n)
		lim_one_side = np.subtract(1, np.power(np.subtract(1, sigm, dtype='float64'), np.divide(1, n + 1, dtype='float64'), dtype='float64'), dtype='float64')
		initialValues['lower'] = a0 if a0 != 0 else 0 if k == 0 else max(lim_binomial, 1.0 / float(n)) + 1e-6
		initialValues['upper'] = b0 if b0 != 0 else 0 if k == n else max(lim_binomial, 1.0 / float(n)) + 1e-6
		limInf['upper'] = 0 if k == n else min(lim_one_side, lim_binomial)
		limInf['lower'] = 0 if k == 0 else min(lim_one_side, lim_binomial)
		limSup['upper'] = 1e-15 if k == n else max(1.0 / float(n), lim_binomial * 2)
		limSup['lower'] = 1e-15 if k == 0 else max(1.0 / float(n), lim_binomial * 2)
		# limSup['upper'] = 0 if k == n else 1.0 - float(k) / float(n)
		# limSup['upper'] = 0 if k == n else 1.0 - float(k) / float(n)
		amin, edm, errdef = ro.Double(0.), ro.Double(0.), ro.Double(0.)
		nvpar, nparx, icstat = ro.Long(0), ro.Long(0), ro.Long(0)
		p, pe = ro.Double(0), ro.Double(0)
		if k != n and k!=0:
			myMinuit.mnparm(0, 'lower', float(initialValues['lower']), 1e-6, float(limInf['lower']), float(limSup['lower']), ierflg)
			if ierflg != 0:
				print 'There was an error starting the lower parameter! ierflg =', ierflg
			if k == 0:
				myMinuit.FixParameter(0)
			myMinuit.mnparm(1, 'upper', float(initialValues['upper']), 1e-6, float(limInf['upper']), float(limSup['upper']), ierflg)
			if ierflg != 0:
				print 'There was an error starting the upper parameter! ierflg =', ierflg
			if k == n:
				myMinuit.FixParameter(1)
			myMinuit.mnparm(2, 'lambd', 0, 1e-6, 0, 0, ierflg)
			if ierflg != 0:
				print 'There was an error starting the lagrange multiplier parameter! ierflg =', ierflg
			myMinuit.mnparm(3, 'k', k, k / float(10) + 0.001, 0, 0, ierflg)
			if ierflg != 0:
				print 'There was an error starting the fixed parameter k! ierflg =', ierflg
			myMinuit.FixParameter(3)
			myMinuit.mnparm(4, 'N', n, n / float(10) + 0.001, 0, 0, ierflg)
			if ierflg != 0:
				print 'There was an error starting the fixed parameter N! ierflg =', ierflg
			myMinuit.FixParameter(4)
			myMinuit.mnparm(5, 'sigma', sigm, sigm / float(10), 0, 0, ierflg)
			if ierflg != 0:
				print 'There was an error starting the fixed parameter sigma! ierflg =', ierflg
			myMinuit.FixParameter(5)
			# Configure Minuit
			arglist = array('d', (0, 0))
			arglist[0] = 2
			myMinuit.mnexcm("SET STRATEGY", arglist, 1, ierflg)
			arglist[0] = max_iter
			arglist[1] = tolerance
			# myMinuit.mnexcm("MINIMIZE", arglist, 2, ierflg)
			# myMinuit.mnexcm("SIMPLEX", arglist, 2, ierflg)
			myMinuit.mnexcm(minimizer, arglist, 2, ierflg)
			# myMinuit.mnexcm("IMPROVE", arglist, 1, ierflg)
			if ierflg != 0:
				print 'There was an error configuring the minimizer! ierflg =', ierflg
			# Initialize Minuit status variables
			# Run Minuit
			myMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat)
			# Get results
			myMinuit.GetParameter(0, p, pe)
			low, lowerr = deepcopy(p), deepcopy(pe)
			myMinuit.GetParameter(1, p, pe)
			up, uperr = deepcopy(p), deepcopy(pe)
			myMinuit.GetParameter(2, p, pe)
			lambd, lambderr = deepcopy(p), deepcopy(pe)
		else:
			low = 0 if k == 0 else np.subtract(1, np.power(np.subtract(1, sigm, dtype='float64'), np.divide(1, n + 1, dtype='float64'), dtype='float64'), dtype='float64')
			up = np.subtract(1, np.power(np.subtract(1, sigm, dtype='float64'), np.divide(1, n + 1, dtype='float64'), dtype='float64'), dtype='float64') if k == 0 else 0
			lambd, lowerr, uperr = 0, 0, 0
		if is_last:
			alt = np.sqrt(k * (1 - float(k) / float(n))) / float(n) + 1e-6
			# low = 0 if k == 0 else alt if low == 0 else low
			low = 0 if k == 0 else low
			# up = 0 if k == n else alt if up == 0 else up
			up = 0 if k == n else up
			return {'lower': low, 'upper': up, 'lambda': lambd, 'lower_err': lowerr, 'upper_err': uperr}
		else:
			return self.FindAsymmetricUncertaintiesWithMinuit(k, n, sigm, max_iter, tolerance, low, up, is_last=True)

	def FindAsymmetricUncertaintiesWithDiscrete(self, k, n, sigm):
		"""
		Numerical method to calculate the asymmetric uncertainties. Is the method used
		:param k: number of successes
		:param n: total number of events
		:param sigm: error bar used for single sided cases
		:return: a dictionary containing the lower error bar, and the upper error bar calculaed using the beta distribution function from scipy stats packate (more precise than root and faster)
		"""
		if k == 0 or k == n:
			edge_value = np.subtract(1, np.power(np.subtract(1, sigm, dtype='float64'), np.divide(1, n + 1, dtype='float64'), dtype='float64'), dtype='float64')
			low = 0 if k == 0 else edge_value
			up = 0 if k == n else edge_value
		else:
			subdiv = self.efficiency_subdiv
			eps_vect = np.append(np.arange(0, 1, np.divide(1, subdiv * n, dtype='float64'), dtype='float64'), 1)
			beta_pdf = sps.beta(k + 1, n - k + 1)
			prob_vect = beta_pdf.pdf(eps_vect).astype('float64')
			# prob_vect = np.array([ro.TMath.BetaDist(eps, k + 1, n - k + 1) for it, eps in np.ndenumerate(eps_vect)], dtype='float64')
			biggest_bins_pos = prob_vect.argsort()[::-1]
			prob_vect_ordered = np.sort(prob_vect)[::-1]
			ordered_integral = np.array([np.divide(prob_vect_ordered[:it].sum(), n * subdiv, dtype='float64') for it in xrange(1, len(prob_vect_ordered) + 1)], dtype='float64')
			num_bins_used = (ordered_integral - sigm >= 0).argmax()
			# ipdb.set_trace()
			eps_in_integral = eps_vect[biggest_bins_pos[:num_bins_used + 1]]
			low, up = np.subtract(eps_in_integral[0], eps_in_integral.min(), dtype='float64'), np.subtract(eps_in_integral.max(), eps_in_integral[0], dtype='float64')
		return {'lower': low, 'upper': up}

	def FitLanGaus(self, name, conv_steps=100, color=ro.kRed, xmin=-10000000, xmax=-10000000):
		"""
		Method used to fit a langaus in a given histogram
		:param name: name of the histogram to be fitted
		:param conv_steps: number of convolution steps for convolving the landau with the gaussian
		:param color: color of the fitted graph
		:param xmin: minimum value for the fit. if -10000000, it will use the limits of the histogram
		:param xmax: maximum value for the fit. if -10000000, it will use the limits of the histogram
		:return:
		"""
		self.canvas[name].cd()
		self.langaus[name] = LanGaus(self.histo[name])
		self.langaus[name].LanGausFit(conv_steps, xmin=xmin, xmax=xmax)
		xlow, xhigh = self.langaus[name].fit_range['min'], self.langaus[name].fit_range['max']
		self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
		self.line[name].SetLineColor(ro.kViolet + 1)
		self.line[name].SetLineWidth(4)
		fitmean = self.langaus[name].fit.Mean(xlow, xhigh)
		self.langaus[name].fit.Draw('same')
		self.langaus[name].fit.SetLineColor(color)
		self.line[name].Draw('same')
		self.histo[name].FindObject('stats').SetOptFit(1)
		self.histo[name].FindObject('stats').SetX1NDC(0.6)
		self.histo[name].FindObject('stats').SetX2NDC(0.9)
		self.histo[name].FindObject('stats').SetY1NDC(0.6)
		self.histo[name].FindObject('stats').SetY2NDC(0.9)
		ro.gPad.Update()
		ps = AddLineToStats(self.canvas[name], 'stats', 'mystats', 'Mean_{Fit}', fitmean)
		self.histo[name].SetStats(0)
		self.canvas[name].Modified()
		ps.Draw()
		self.trash.append(ps)
		ro.gPad.Update()
	# print '{n}: <PH> ex= {f}'.format(n=name, f=fitmean)

	def DrawDoubleLangaus(self, name, name1, name2, color=ro.kBlack):
		"""
		Method used to fit the sum of two langaus
		:param name:
		:param name1:
		:param name2:
		:param color:
		:return:
		"""
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
			AddLineToStats(self.canvas[name], 'stats', 'mystats', 'Mean_{Fit}', fitmean)
			self.histo[name].SetStats(0)
			self.canvas[name].Modified()
			ro.gPad.Update()
			print '{n}: <PH> = {f}'.format(n=name, f=fitmean)

	def TwoLanGaus(self, x, params):
		"""
		Function for ROOT of the addition of two langaus
		:param x:
		:param params:
		:return:
		"""
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

	def FitGaus(self, name, num_sigma=2, iterations=2):
		"""
		Method used to fit a gaussian on a desired histogram. The fit, will be iterated taking only the limits specified by num_sigma within the mean value of the fit
		:param name: name of the histogram to be fitted
		:param num_sigma: number of sigmas away from the mean value used for the fit.
		:param iterations: number of iterations to do the fit. Each iteration, the limits of the fit change
		:return:
		"""
		if name in self.histo.keys():
			if name in self.canvas.keys():
				self.canvas[name].cd()
			# histo = self.histo[name]
			xmean, xrms = self.histo[name].GetMean(), self.histo[name].GetRMS()
			xmin, xmax = xmean - num_sigma * xrms, xmean + num_sigma * xrms
			func = ro.TF1('f_gaus_' + name, 'gaus', xmean - (num_sigma + 1) * xrms, xmean + (num_sigma + 1) * xrms)
			func.SetNpx(1000)
			func.SetLineStyle(1)
			func.SetLineColor(ro.kRed)
			func.SetLineWidth(2)
			params = np.array((self.histo[name].GetBinContent(self.histo[name].GetXaxis().FindBin(xmean)), xmean, xrms), 'float64')
			func.SetParameters(params)
			for it in xrange(iterations):
				temp = self.histo[name].Fit('f_gaus_' + name, 'QIEBMSN', 'goff', xmin, xmax)
				params = np.array((temp.Parameter(0), temp.Parameter(1), temp.Parameter(2)), 'float64')
				func.SetParameters(params)
			self.fits[name] = self.histo[name].Fit('f_gaus_' + name, 'QIEBMS', 'goff', params[1] - num_sigma * params[2], params[1] + num_sigma * params[2])
			self.histo[name].GetFunction('f_gaus_' + name).Draw('same')
			ro.gPad.Update()

	def FitPol(self, name, num_pol=0):
		"""
		Method used to fit a polynomial function of degree num_pol
		:param name: name of the histogram to be fitted
		:param num_pol: degree of the polynomial to be fitted
		:return:
		"""
		if name in self.histo.keys():
			if name in self.canvas.keys():
				self.canvas[name].cd()
			xmin, xmax = self.histo[name].GetXaxis().GetBinLowEdge(self.histo[name].FindFirstBinAbove()), self.histo[name].GetXaxis().GetBinLowEdge(self.histo[name].FindLastBinAbove() + 1)
			func = ro.TF1('f_pol{num}_{n}'.format(num=num_pol, n=name), 'pol{n}'.format(n=num_pol), xmin, xmax)
			func.SetNpx(1000)
			func.SetLineStyle(1)
			func.SetLineColor(ro.kRed)
			func.SetLineWidth(2)
			self.fits[name] = self.histo[name].Fit('f_pol{num}_{n}'.format(num=num_pol, n=name), 'QIEBMS', 'goff', xmin, xmax)
			self.histo[name].GetFunction('f_pol{num}_{n}'.format(num=num_pol, n=name)).Draw('same')
			ro.gPad.Update()

	def SaveCanvasInlist(self, list):
		"""
		Method used to save the plots with names in the given list
		:param list: list of names of the plots to be saved
		:return:
		"""
		if not os.path.isdir('{d}/{r}/{sd}'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/{sd}'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir))
		for canvas in list:
			if self.canvas.has_key(canvas):
				if self.canvas[canvas]:
					self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.png'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))
					# self.canvas[canvas].Print('{d}/{r}/{sd}/{c}.png'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))
					self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.root'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))
					# self.canvas[canvas].Print('{d}/{r}/{sd}/{c}.root'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))

	def LoadPlotsInSubdir(self):
		"""
		Try to load all the root files of the plots previously saved
		:return:
		"""
		if not os.path.isdir('{d}/{r}/{sd}'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir)):
			print 'The working subdirectory: {d}/{r}/{sd} does not exist'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir)
			return
		root_files_list = glob('{d}/{r}/{sd}/*.root'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir))
		for file_path_i in root_files_list:
			file_i = ro.TFile(file_path_i, 'read')
			canvas_name = file_i.GetListOfKeys()[0].GetTitle()
			family_name = canvas_name[2:] if canvas_name.startswith('c_') else canvas_name
			canvas = file_i.Get(canvas_name)
			self.canvas[family_name] = canvas
			histo = canvas.GetPrimitive('h_' + family_name)
			if histo:
				if self.IsObjectTH(histo):
					self.histo[family_name] = histo
				if self.IsObjectTProfile2D(histo):
					self.profile[family_name] = histo
			file_i.Close()

	def IsObjectTH(self, obj):
		"""
		Check if the object is of type TH1 or TH2
		:param obj: object to be checked
		:return: returns True if the object is of type TH1 or TH2
		"""
		if obj:
			if obj.InheritsFrom(ro.TH1.Class().GetName()) or obj.InheritsFrom(ro.TH2.Class().GetName()):
				return True
		return False

	def IsObjectTProfile2D(self, obj):
		"""
		Check if the object is of type TProifile2D
		:param obj: object to be checked
		:return: returns True if the object is of type TProfile2D
		"""
		if obj:
			if obj.InheritsFrom(ro.TProfile2D.Class().GetName()):
				return True
		return False

	def CreateFriendWithSaturationRegions(self, suffix='', skipAfter=0, skipBefore=0):
		"""
		Creates a tree friend with saturation regions flag specified by skipAfter and skipBefore.
		The normal case, is to flag only the saturated event, which corresponds to skipAfter = 1, skipBefore = 0
		:param suffix: suffix used for saving the tree. This enables the possibility to have multiple files in the sub directory
		:param skipAfter: number of events to skip after the saturation event including the saturated event
		:param skipBefore: number of events to skip before the saturation event without including the saturated event
		:return:
		"""
		if not self.cuts_man.sat_evts:
			self.cuts_man.FindSaturationEvents()
		satFile = ro.TFile('{d}/{r}/satRegions{s}.{r}.root'.format(d=self.dir, s=suffix, r=self.run), 'RECREATE')
		satTree = ro.TTree('satRegions', 'satRegions')
		satRegEv = np.zeros(1, '?')
		satTree.Branch('satRegion', satRegEv, 'satRegion/O')
		ev0 = int(self.trans_tree.GetMinimum('event'))
		evMax = int(self.trans_tree.GetMaximum('event') + 1)
		nevents = evMax - ev0
		bar = CreateProgressBarUtils(nevents)
		bar.start()
		for ev in xrange(ev0, evMax):
			satRegEv.fill(0)
			if self.cuts_man.IsEventInSaturationRegion(ev, skipAfter=skipAfter, skipBefore=skipBefore):
				satRegEv.fill(1)
			satTree.Fill()
			bar.update(ev - ev0 + 1)
		satFile.Write()
		satFile.Close()
		bar.finish()

	def AddFriendWithSaturationRegions(self, skipAfter=100, skipBefore=0):
		"""
		Adds a friend tree with the saturation flag corresponding to the specified number of events before and after saturation.
		It will call CreateFriendWithSaturationRegions if the tree does not exist
		:param skipAfter: number of events to skip after the saturation event including the saturated event
		:param skipBefore: number of events to skip before the saturation event without including the saturated event
		:return:
		"""
		suffix = '{sb}b{sa}a'.format(sb=skipBefore, sa=skipAfter)
		if not self.trans_tree.GetFriend('satRegions'):
			if os.path.isfile('{d}/{r}/satRegions{s}.{r}.root'.format(d=self.dir, s=suffix, r=self.run)):
				self.trans_tree.AddFriend('satRegions', '{d}/{r}/satRegions{s}.{r}.root'.format(d=self.dir, s=suffix, r=self.run))
			else:
				self.CreateFriendWithSaturationRegions(suffix, skipAfter, skipBefore)
				self.AddFriendWithSaturationRegions(skipAfter=skipAfter, skipBefore=skipBefore)

	def CreateFriendWithNewPedestalBuffer(self, slide_length=50, hit_factor=0, seed_factor=0):
		"""
		Creates a new tree friend with a different pedestal calculation due to a different buffer length or hit, seed factors
		:param slide_length: buffer length used for the calculation. The original buffer in the analysis is 500. It is recommended to use values >=~ 200
		:param hit_factor: hit factor to calculate the pedestal.
		:param seed_factor: seed factor used to calculate the pedestal
		:return:
		"""
		hit_fact = hit_factor if hit_factor != 0 else self.hit_factor
		seed_fact = seed_factor if seed_factor != 0 else self.seed_factor
		if hit_fact + seed_fact == 0:
			ExitMessage('Hit and seed factor are 0. Exiting', os.EX_NOINPUT)
		ev_ini, ev_end = self.trans_tree.GetMinimum('event'), self.trans_tree.GetMaximum('event')
		self.OpenPedestalFileAndTree()
		pedCalc = PedestalCalculations(self.ped_tree, self.dir, self.run, slide_length, hit_fact, seed_fact, ev_ini, ev_end)
		self.CloseOriginalPedestalFile()
		pedCalc.CalculateDevicesPedestals()

	def AddFriendWithNewPedestalBuffer(self, slide_length=50, hit_factor=0, seed_factor=0):
		"""
		Tries to add a tree friend with a different calculation of the pedestal. If the tree does not exist, the method will create it by calling CreateFriendWithNewPedestalBuffer
		:param slide_length: buffer length used for the calculation. The original buffer in the analysis is 500. It is recommended to use values >=~ 200
		:param hit_factor: hit factor to calculate the pedestal.
		:param seed_factor: seed factor used to calculate the pedestal
		:return:
		"""
		hit_fact = hit_factor if hit_factor != 0 else self.hit_factor
		seed_fact = seed_factor if seed_factor != 0 else self.seed_factor
		if not self.trans_tree.GetFriend('pedTree'):
			if os.path.isfile('{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.dir, s=slide_length, r=self.run)):
				self.trans_tree.AddFriend('pedTree', '{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.dir, r=self.run, s=slide_length))
				self.noise_friend_buffer = slide_length
			else:
				self.CreateFriendWithNewPedestalBuffer(slide_length, hit_fact, seed_fact)
				if os.path.isfile('{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.dir, s=slide_length, r=self.run)):
					self.trans_tree.AddFriend('pedTree', '{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.dir, r=self.run, s=slide_length))
					self.noise_friend_buffer = slide_length
				else:
					print 'Something went wrong... Created file does not exist!'

	def UnfriendTree(self, extreefriend):
		"""
		Method used to unfriend a friend from the transparent tree
		:param extreefriend: is the tree object to be unfriend. It can be given by trans_tree.GetFriend(<friend_name>)
		:return:
		"""
		treename = extreefriend.GetName()
		hasfriend = True if self.trans_tree.GetFriend(treename) else False
		if hasfriend:
			self.trans_tree.RemoveFriend(extreefriend)

	def CheckIfPedTreeFriendExists(self, buff=50):
		"""
		Checks if the pedestal tree friend exists for a certain buffer size
		:param buff: buffer size of the pedestal tree friend
		:return: returns True if it exist. Otherwise False
		"""
		return True if os.path.isfile('{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.dir, s=buff, r=self.run)) else False

	def CreatePedTreeFriendsForStudy(self, buffers_array):
		"""
		Creates several pedestal tree friends in parallel. The number of parallel calculations is given by num_parallel which can be changed before calling this method
		:param buffers_array: numpy array containing the sizes of the buffers of the tree friends to be created
		:return:
		"""
		ev_ini, ev_end = self.trans_tree.GetMinimum('event'), self.trans_tree.GetMaximum('event')
		job_chunks = [buffers_array[i:i+self.num_parallel] for i in xrange(0, buffers_array.size, self.num_parallel)]
		for buffsrow in job_chunks:
			process = []
			for it, buff in enumerate(buffsrow):
				if it == len(buffsrow) - 1:
					self.OpenPedestalFileAndTree()
					tempPedCalc = PedestalCalculations(self.ped_tree, self.dir, self.run, buff, self.hit_factor, self.seed_factor, ev_ini, ev_end, True)
					self.CloseOriginalPedestalFile()
				else:
					self.OpenPedestalFileAndTree()
					tempPedCalc = PedestalCalculations(self.ped_tree, self.dir, self.run, buff, self.hit_factor, self.seed_factor, ev_ini, ev_end, False)
					self.CloseOriginalPedestalFile()
				tempPedCalc.start()
				process.append(tempPedCalc)
			for j in process:
				j.join()
		print ('Finished creating the pedTrees for buffers:', buffers_array) if len(buffers_array) > 0 else ('All pedTrees already exist ')

	def AddFriendWithCells(self, xoff=0, yoff=0, overWrite=False):
		"""
		Tries to add a friend with cell information with a xoff in width_units and a yoff in um. If it is not in the root file, it will try to create one. The tree name has the offset in nm with 100s of nm in resolution
		:param xoff: offset in x in cell's width units
		:param yoff: offset in y in um
		:param overWrite: if True, deletes tree and creates new one
		:return:
		"""
		xoff_nm = int(RoundInt(xoff * self.row_cell_info_diamond['width'] * 10) * 100)
		yoff_nm = int(RoundInt(yoff * 10) * 100)
		x_prefix = '' if xoff >= 0 else 'Neg'
		y_prefix = '' if yoff >= 0 else 'Neg'
		treename = 'cells{xp}{x}X{yp}{y}Y'.format(x=abs(xoff_nm), xp=x_prefix, yp=y_prefix, y=abs(yoff_nm))
		if not self.trans_tree.GetFriend(treename):
			if os.path.isfile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run)):
				if overWrite and xoff==0 and yoff==0:
					os.remove('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run))
					self.ResetAreas()
					self.CreateGridAreas()
					self.CreateTCutGs()
					self.CreateFriendWithCells(xoff, yoff)
					self.AddFriendWithCells(xoff, yoff)
				else:
					tempf = ro.TFile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run), 'READ')
					if tempf.Get(treename):
						tempf.Close()
						self.trans_tree.AddFriend(treename, '{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run))
					else:
						tempf.Close()
						self.CreateFriendWithCells(xoff, yoff)
						self.AddFriendWithCells(xoff, yoff)
			else:
				self.CreateFriendWithCells(xoff, yoff)
				self.AddFriendWithCells(xoff, yoff)
		else:
			if overWrite and xoff==0 and yoff==0:
				self.UnfriendTree(self.trans_tree.GetFriend(treename))
				if os.path.isfile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run)):
					os.remove('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run))
				self.ResetAreas()
				self.CreateGridAreas()
				self.CreateTCutGs()
				self.CreateFriendWithCells(xoff, yoff)
				self.AddFriendWithCells(xoff, yoff)


	def CreateFriendWithCells(self, xoff=0, yoff=0):
		"""
		Creates a tree friend with a specified offset xoff, yoff. The tree will be save in the same Root File with a different tree name
		:param xoff: offset from the center position given by the alignment. Normally is a small number smaller than 0.01. The number is in units of the pitch of the column
		:param yoff: offset from the pickle values '0_even' and '0_odd' to create the TCutGs of the cells. It is given in mum
		:return:
		"""
		if self.dia_cols:
			print 'Getting vectors to calculate cells...', ; sys.stdout.flush()
			tempc = self.cuts_man.transp_ev
			leng = self.trans_tree.Draw('event:diaChXPred:diaChYPred', tempc, 'goff')
			if leng == -1:
				ExitMessage('Error, could not get branches diaChXPred or diaChYPred. Check tree', os.EX_DATAERR)
			while leng > self.trans_tree.GetEstimate():
				self.trans_tree.SetEstimate(leng)
				leng = self.trans_tree.Draw('event:diaChXPred:diaChYPred', tempc, 'goff')
			tempev = self.trans_tree.GetVal(0)
			tempx = self.trans_tree.GetVal(1)
			tempy = self.trans_tree.GetVal(2)
			ev_vect = [int(tempev[ev]) for ev in xrange(leng)]
			x_vect = [tempx[ev] for ev in xrange(leng)]
			y_vect = [tempy[ev] for ev in xrange(leng)]
			print 'Done'
			col_dic = {}
			row_dic = {}
			x0_dic = {}
			y0_dic = {}
			print 'Calculating cells for each transparent event:'
			bar = CreateProgressBarUtils(leng)
			bar.start()
			for ev in xrange(leng):
				dist = np.array([[cell.GetDistanceToCenter(x_vect[ev], y_vect[ev]) for cell in col.cells] for col in self.dia_cols.cols]).flatten()
				cols = np.array([[cell.col_num for cell in col.cells] for col in self.dia_cols.cols]).flatten()
				rows = np.array([[cell.row_num for cell in col.cells] for col in self.dia_cols.cols]).flatten()

				argmin = dist.argmin()
				colmin = cols[argmin]
				rowmin = rows[argmin]

				col_dic[ev_vect[ev]] = colmin
				row_dic[ev_vect[ev]] = rowmin
				x0_dic[ev_vect[ev]] = x_vect[ev] - self.dia_cols.cols[colmin].cells[rowmin].xcenter
				y0_dic[ev_vect[ev]] = y_vect[ev] - self.dia_cols.cols[colmin].cells[rowmin].ycenter
				bar.update(ev + 1)
			bar.finish()

			print 'Saving tree:'
			xoff_nm = int(RoundInt(xoff * self.row_cell_info_diamond['width'] * 10) * 100)
			yoff_nm = int(RoundInt(yoff * 10) * 100)
			x_prefix = '' if xoff >= 0 else 'Neg'
			y_prefix = '' if yoff >= 0 else 'Neg'
			treename = 'cells{xp}{x}X{yp}{y}Y'.format(x=abs(xoff_nm), xp=x_prefix, yp=y_prefix, y=abs(yoff_nm))
			cellsFile = ro.TFile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run), 'UPDATE') if os.path.isfile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run)) else ro.TFile('{d}/{r}/cells.{r}.root'.format(d=self.dir, r=self.run), 'RECREATE')
			cellsTree = ro.TTree(treename, treename)
			x0 = np.zeros(1, 'f4')
			y0 = np.zeros(1, 'f4')
			col = np.zeros(1, 'uint8')
			row = np.zeros(1, 'uint8')
			cellsTree.Branch('x0', x0, 'x0/F')
			cellsTree.Branch('y0', y0, 'y0/F')
			cellsTree.Branch('col', col, 'col/b')
			cellsTree.Branch('row', row, 'row/b')
			ev0 = int(self.trans_tree.GetMinimum('event'))
			evMax = int(self.trans_tree.GetMaximum('event') + 1)
			nevents = evMax - ev0
			bar = CreateProgressBarUtils(nevents)
			bar.start()
			for ev in xrange(ev0, evMax):
				if ev in ev_vect:
					x0.fill(x0_dic[ev])
					y0.fill(y0_dic[ev])
					col.fill(col_dic[ev])
					row.fill(row_dic[ev])
				else:
					x0.fill(0)
					y0.fill(0)
					col.fill(0)
					row.fill(0)
				cellsTree.Fill()
				bar.update(ev - ev0 + 1)
			cellsFile.Write()
			cellsFile.Close()
			bar.finish()
		else:
			ExitMessage('Cant create friend with cell info. Check', os.EX_SOFTWARE)

	def CloseInputROOTFiles(self):
		if self.trans_file:
			if self.trans_file.IsOpen():
				self.trans_file.Close()
			if self.trans_tree:
				del self.trans_tree
			# self.trans_file = None
			# self.trans_tree = None

	def CloseOriginalPedestalFile(self):
		if self.ped_file:
			if self.ped_file.IsOpen():
				self.ped_file.Close()
			if self.ped_tree:
				del self.ped_tree


	def GetPHChVar(self, ch, typ='H', isSNR=False, isFriend=False):
		"""
		Method that returns the PH variable used for a certain number of channel
		:param ch: The number of the channel considered for the calculation
		:param typ: Type of the PH variable. If 'H' will take into account the highest channels. If 'Ch', it will take into account the channels ordered using the proximity from the predicted hit position as criteria.
		:param isSNR: flag that enables the option to obtain the variable in terms of Sigma (SNR)
		:param isFriend: Flag used to indicate that the signal, or noise values will be used from a previously befriended pedestal tree (pedTree) instead to the original transparent Tree values.
		:return: returns the string that is used in the analysis as a PH variable
		"""
		varch = 'clusterChannel' + str(ch) if typ == 'Ch' else 'clusterChannelHighest' + str(ch) if typ == 'H' else str(ch)
		return '(pedTree.diaChSignal[{c}])'.format(c=varch) if isFriend and not isSNR else '(diaChSignal[{c}]/diaChPedSigmaCmc[{c}])'.format(c=varch) if not isFriend and isSNR else '(pedTree.diaChSignal[{c}]/pedTree.diaChPedSigmaCmc[{c}])'.format(c=varch) if isFriend and isSNR else '(diaChSignal[{c}])'.format(c=varch)

	def GetPHNChsVar(self, n, typ='H', isSNR=False, isFriend=False):
		"""
		Method that returns the cumulative PH variable used for a certain number of channels
		:param ch: The number of the cumulative channels considered for the calculation
		:param typ: Type of the PH variable. If 'H' will take into account the highest channels. If 'Ch', it will take into account the channels ordered using the proximity from the predicted hit position as criteria.
		:param isSNR: flag that enables the option to obtain the variable in terms of Sigma (SNR)
		:param isFriend: Flag used to indicate that the signal, or noise values will be used from a previously befriended pedestal tree (pedTree) instead to the original transparent Tree values.
		:return: returns the string that is used in the analysis as a PH variable
		"""
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetPHChVar(chi - 1 if typ == 'Ch' else chi, typ, isSNR, isFriend))
		return '+'.join(temp)

	def FindMaxMinVarz(self):
		"""
		min max variables used for each plotting case... should be removed as phmin and phmax is more universal and this lags the analysis. The calculation is saved in a pickle to prevent future recalculations
		:return:
		"""
		if os.path.isfile('{d}/{r}/{s}/min_max_values.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
			with open('{d}/{r}/{s}/min_max_values.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir), 'rb') as pkl:
				minmax_temp = pickle.load(pkl)
				self.minz, self.maxz = minmax_temp['min'], minmax_temp['max']
				return
		print 'Finding Maximum and Minimum ranges for plotting...', ; sys.stdout.flush()
		for ch in xrange(self.cluster_size):
			self.minz['PH_Ch{i}_adc'.format(i=ch)] = GetMinimumFromTree(self.trans_tree, self.GetPHChVar(ch, 'Ch'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_Ch{i}'.format(i=ch))]))
			self.maxz['PH_Ch{i}_adc'.format(i=ch)] = GetMaximumFromTree(self.trans_tree, self.GetPHChVar(ch, 'Ch'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_Ch{i}'.format(i=ch))]))
			self.minz['PH_H{i}_adc'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHChVar(ch+1, 'H'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_H{i}'.format(i=ch+1))]))
			self.maxz['PH_H{i}_adc'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHChVar(ch+1, 'H'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_H{i}'.format(i=ch+1))]))
			self.minz['PH_Ch{i}_snr'.format(i=ch)] = GetMinimumFromTree(self.trans_tree, self.GetPHChVar(ch, 'Ch', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_Ch{i}'.format(i=ch))]))
			self.maxz['PH_Ch{i}_snr'.format(i=ch)] = GetMaximumFromTree(self.trans_tree, self.GetPHChVar(ch, 'Ch', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_Ch{i}'.format(i=ch))]))
			self.minz['PH_H{i}_snr'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHChVar(ch+1, 'H', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_H{i}'.format(i=ch+1))]))
			self.maxz['PH_H{i}_snr'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHChVar(ch+1, 'H', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH_H{i}'.format(i=ch+1))]))
			self.minz['PH{i}_Ch_adc'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'Ch'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_Ch'.format(i=ch+1))]))
			self.maxz['PH{i}_Ch_adc'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'Ch'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_Ch'.format(i=ch+1))]))
			self.minz['PH{i}_Ch_snr'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'Ch', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_Ch'.format(i=ch+1))]))
			self.maxz['PH{i}_Ch_snr'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'Ch', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_Ch'.format(i=ch+1))]))
			self.minz['PH{i}_H_adc'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'H'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_H'.format(i=ch+1))]))
			self.maxz['PH{i}_H_adc'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'H'), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_H'.format(i=ch+1))]))
			self.minz['PH{i}_H_snr'.format(i=ch+1)] = GetMinimumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'H', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_H'.format(i=ch+1))]))
			self.maxz['PH{i}_H_snr'.format(i=ch+1)] = GetMaximumFromTree(self.trans_tree, self.GetPHNChsVar(ch+1, 'H', True), self.cuts_man.AndCuts([self.cuts_man.transp_ev, self.cuts_man.GetPHCuts('PH{i}_H'.format(i=ch+1))]))
		print 'Done'
		minmax_obj = {'min': self.minz, 'max': self.maxz}
		if not os.path.isdir('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir)):
			os.makedirs('{d}/{r}/{s}'.format(d=self.dir, r=self.run, s=self.pkl_sbdir))
		print 'Saving limits in pickle file', '{d}/{r}/{s}/min_max_values.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir), '...', ; sys.stdout.flush()
		pickle.dump(minmax_obj, open('{d}/{r}/{s}/min_max_values.{r}.pkl'.format(d=self.dir, r=self.run, s=self.pkl_sbdir), 'wb'))
		print 'Done'

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
	parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
	parser.add_option('-c', '--colpitch', dest='colpitch', type='int', default=50, help='column pitch of the device')
	parser.add_option('-n', '--numstrips', dest='numstrips', type='int', default=2, help='Number of strips to use')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	numstrips = int(options.numstrips)
	colpitch = int(options.colpitch)

	tg = TransparentGrid(dir=dir, run=run, col_pitch=colpitch)
