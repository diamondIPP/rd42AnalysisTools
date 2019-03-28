#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class CenterCellAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize):
		self.phdelta = 0
		self.window_shift = 3
		# self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.125
		self.delta_ev = 100
		self.min_snr, self.max_snr = -650, 650
		self.min_adc, self.max_adc = -6500, 6500
		self.delta_adc, self.delta_snr = 20, 2
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.trash = []
		self.w = 0
		self.trans_grid = trans_grid
		self.num_strips = numstrips
		self.cluster_size = clustersize

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.in_transp_cluster = ''

		self.noise_varz = {}
		self.noise_friend_varz = {}

		self.ph_adc_ch_varz = {}
		self.ph_snr_ch_varz = {}
		self.ph_adc_h_varz = {}
		self.ph_snr_h_varz = {}
		self.phN_adc_ch_varz = {}
		self.phN_snr_ch_varz = {}
		self.phN_adc_h_varz = {}
		self.phN_snr_h_varz = {}

		self.ph_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}

		self.sat_adc_ch_cut = {}
		self.sat_adc_h_cut = {}
		self.not_sat_adc_ch_cut = {}
		self.not_sat_adc_h_cut = {}

		self.sat_adc_N_ch_cut = {}
		self.sat_adc_N_h_cut = {}
		self.not_sat_adc_N_ch_cut = {}
		self.not_sat_adc_N_h_cut = {}


		self.in_central_rect_cut = {}
		self.out_central_rect_cut = {}

		self.sat_adc_inside_cut = {}
		self.sat_adc_outside_cut = {}
		self.nosat_adc_inside_cut = {}
		self.nosat_adc_outside_cut = {}

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def DoCenterRegionStudies(self, p_array, cells='all'):
		self.GetCutsFromCutManager(cells)
		self.GetVarzFromTranspGrid()
		self.DoPlots(p_array, cells)

	def DoPlots(self, p_array, cells='all'):
		suffix = self.suffix[cells]
		for p in p_array:
			if 0 < p < 100:
				varzkey = 'PH2_H'
				if p in self.in_central_rect_cut.keys() and p in self.out_central_rect_cut.keys():
					tempcutin = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_h_cuts[cells][varzkey], self.in_central_rect_cut[p])
					tempcutout = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_h_cuts[cells][varzkey], self.out_central_rect_cut[p])

					self.DoCellOverlayPlot('PH2_H_adc_cells_in_rect_{p}_percent_{s}'.format(p=p, s=suffix), self.phN_adc_h_varz[varzkey], 'PH2 highest chs [ADC]', tempcutin, cells)
					self.trans_grid.DrawCentralArea('PH2_H_adc_cells_in_rect_{p}_percent_{s}'.format(p=p, s=suffix), p)
					self.DoPHPlot('PH2_H_adc_in_rect_{p}_percent_{s}'.format(p=p, s=suffix),self.phN_adc_h_varz[varzkey], 'PH2 highest chs [ADC]', tempcutin, cells)

					self.DoCellOverlayPlot('PH2_H_adc_cells_out_rect_{p}_percent_{s}'.format(p=p, s=suffix), self.phN_adc_h_varz[varzkey], 'PH2 highest chs [ADC]', tempcutout, cells)
					self.trans_grid.DrawCentralArea('PH2_H_adc_cells_out_rect_{p}_percent_{s}'.format(p=p, s=suffix), p)
					self.DoPHPlot('PH2_H_adc_out_rect_{p}_percent_{s}'.format(p=p, s=suffix),self.phN_adc_h_varz[varzkey], 'PH2 highest chs [ADC]', tempcutout, cells)

	def DoCellOverlayPlot(self, name, varz, varname, cut, cells='all'):
		self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname)
		self.PosCanvas(name)

	def DoPHPlot(self, name, varz, varname, cut, cells='all'):
		self.trans_grid.DrawPHInArea(name, varz, cells, cut, varname=varname)
		self.PosCanvas(name)

	def GetCutsFromCutManager(self, cells):
		self.in_central_rect_cut = self.trans_grid.cuts_man.in_central_rect_region
		self.out_central_rect_cut = self.trans_grid.cuts_man.out_central_rect_region

		self.ph_adc_ch_cuts[cells] = self.trans_grid.cuts_man.ph_adc_ch_cuts[cells]
		self.ph_snr_ch_cuts[cells] = self.trans_grid.cuts_man.ph_snr_ch_cuts[cells]
		self.ph_adc_h_cuts[cells] = self.trans_grid.cuts_man.ph_adc_h_cuts[cells]
		self.ph_snr_h_cuts[cells] = self.trans_grid.cuts_man.ph_snr_h_cuts[cells]
		self.phN_adc_ch_cuts[cells] = self.trans_grid.cuts_man.phN_adc_ch_cuts[cells]
		self.phN_snr_ch_cuts[cells] = self.trans_grid.cuts_man.phN_snr_ch_cuts[cells]
		self.phN_adc_h_cuts[cells] = self.trans_grid.cuts_man.phN_adc_h_cuts[cells]
		self.phN_snr_h_cuts[cells] = self.trans_grid.cuts_man.phN_snr_h_cuts[cells]

		self.sat_adc_ch_cut = self.trans_grid.cuts_man.sat_adc_ch
		self.sat_adc_h_cut = self.trans_grid.cuts_man.sat_adc_h
		self.not_sat_adc_ch_cut = self.trans_grid.cuts_man.not_sat_adc_ch
		self.not_sat_adc_h_cut = self.trans_grid.cuts_man.not_sat_adc_h
		self.sat_adc_N_ch_cut = self.trans_grid.cuts_man.sat_adc_N_ch
		self.sat_adc_N_h_cut = self.trans_grid.cuts_man.sat_adc_N_h
		self.not_sat_adc_N_ch_cut = self.trans_grid.cuts_man.not_sat_adc_N_ch
		self.not_sat_adc_N_h_cut = self.trans_grid.cuts_man.not_sat_adc_N_h

		self.SetCutsForSaturationRatio(cells)

	def GetVarzFromTranspGrid(self):
		self.noise_varz = self.trans_grid.noise_varz
		self.noise_friend_varz = self.trans_grid.noise_friend_varz
		self.ph_adc_ch_varz = self.trans_grid.ph_adc_ch_varz
		self.ph_snr_ch_varz = self.trans_grid.ph_snr_ch_varz
		self.ph_adc_h_varz = self.trans_grid.ph_adc_h_varz
		self.ph_snr_h_varz = self.trans_grid.ph_snr_h_varz
		self.phN_adc_ch_varz = self.trans_grid.phN_adc_ch_varz
		self.phN_snr_ch_varz = self.trans_grid.phN_snr_ch_varz
		self.phN_adc_h_varz = self.trans_grid.phN_adc_h_varz
		self.phN_snr_h_varz = self.trans_grid.phN_snr_h_varz

	def SetCutsForSaturationRatio(self, cells='all'):
		tempsat = self.trans_grid.cuts_man.ConcatenateCuts(self.trans_grid.cuts_man.transp_ev, self.sat_adc_N_h_cut['2_H'])
		tempsat = self.trans_grid.cuts_man.ConcatenateCuts(tempsat, self.phN_adc_h_cuts[cells]['PH2_H'])
		tempnosat = self.trans_grid.cuts_man.ConcatenateCuts(self.trans_grid.cuts_man.transp_ev, self.not_sat_adc_N_h_cut['2_H'])
		tempnosat = self.trans_grid.cuts_man.ConcatenateCuts(tempnosat, self.phN_adc_h_cuts[cells]['PH2_H'])
		for p, cut in self.in_central_rect_cut.iteritems():
			self.sat_adc_inside_cut[p] = self.trans_grid.cuts_man.ConcatenateCuts(tempsat, cut)
			self.nosat_adc_inside_cut[p] = self.trans_grid.cuts_man.ConcatenateCuts(tempnosat, cut)
		for p, cut in self.out_central_rect_cut.iteritems():
			self.sat_adc_outside_cut[p] = self.trans_grid.cuts_man.ConcatenateCuts(tempsat, cut)
			self.nosat_adc_outside_cut[p] = self.trans_grid.cuts_man.ConcatenateCuts(tempnosat, cut)


if __name__ == '__main__':
	c = CenterCellAnalysis(None, 0, 0)
