#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class FinalAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana=None):
		self.phdelta = 0
		self.window_shift = 3
		self.min_snr_neg, self.max_snr_neg = -64, 0
		self.min_adc_neg, self.max_adc_neg = -650, 0
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
		self.noise_ana = noise_ana

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.ph_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}

		self.minz = {t: {} for t in ['all', 'good', 'bad']}
		self.maxz = {t: {} for t in ['all', 'good', 'bad']}

		self.sat_adc_ch_cut = {}
		self.sat_adc_h_cut = {}
		self.not_sat_adc_ch_cut = {}
		self.not_sat_adc_h_cut = {}
		self.sat_adc_N_ch_cut = {}
		self.sat_adc_N_h_cut = {}
		self.not_sat_adc_N_ch_cut = {}
		self.not_sat_adc_N_h_cut = {}

		self.sat_evts_region = ''
		self.not_sat_evts_region = ''

		self.in_transp_cluster = ''

		self.noise_varz = {}
		self.noise_friend_varz = {}
		
		self.neg_snr_ph_ch = {}
		self.neg_snr_ph_h = {}
		self.not_neg_snr_ph_ch = {}
		self.not_neg_snr_ph_h = {}

		self.neg_adc_ph_ch = {}
		self.neg_adc_ph_h = {}
		self.not_neg_adc_ph_ch = {}
		self.not_neg_adc_ph_h = {}

		self.neg_snr_phN_ch = {}
		self.neg_snr_phN_h = {}
		self.not_neg_snr_phN_ch = {}
		self.not_neg_snr_phN_h = {}

		self.neg_adc_phN_ch = {}
		self.neg_adc_phN_h = {}
		self.not_neg_adc_phN_ch = {}
		self.not_neg_adc_phN_h = {}

		self.noise_varz = {}
		self.ph_adc_h_varz = {}
		self.ph_adc_ch_varz = {}
		self.ph_snr_h_varz = {}
		self.ph_snr_ch_varz = {}

		self.phN_adc_h_varz = {}
		self.phN_adc_ch_varz = {}
		self.phN_snr_h_varz = {}
		self.phN_snr_ch_varz = {}

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def DoDeviceMaps(self, cells, cuts='', suffix='no_cuts'):
		def DrawProfile2D(name, varz, varzname, cut, getOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondMap(name, varz, varzname, cells, cut, 'prof colz')
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)
			if getOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.trans_grid.histo['hit_map_' + name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_' + name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
				self.trans_grid.DrawTCutGs('hit_map_' + name, 'diamond')
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_' + name)
				self.PosCanvas('hit_map_' + name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		for ch in xrange(1, self.cluster_size + 1):
			DrawProfile2D('PH{c}_Ch_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [SNR]'.format(c=ch), tempcsnr, False)
			DrawProfile2D('PH{c}_Ch_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [ADC]'.format(c=ch), tempcadc, False)
			DrawProfile2D('PH{c}_H_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [SNR]'.format(c=ch), tempcsnr, ch == self.cluster_size)
			DrawProfile2D('PH{c}_H_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [ADC]'.format(c=ch), tempcadc, ch == self.cluster_size)

	def DoCellMaps(self, cells, cuts='', suffix='no_cuts'):
		def PlotCellsProfiles(name, varz, zmin, zmax, varname, cut, doOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname)
			self.trans_grid.profile[name].SetMinimum(min(0, zmin))
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.PosCanvas(name)
			if doOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.PosCanvas('hit_map_' + name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		for ch in xrange(1, self.cluster_size + 1):
			minz, maxz = self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)]
			PlotCellsProfiles('PH{c}_Ch_cell_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], minz, max(maxz, 480), 'PH{c} cluster chs [SNR]', tempcsnr)
			minz, maxz = self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)]
			PlotCellsProfiles('PH{c}_Ch_cell_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], minz, max(maxz, 4800), 'PH{c} cluster chs [ADC]', tempcadc)
			minz, maxz = self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)]
			PlotCellsProfiles('PH{c}_H_cell_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], minz, max(maxz, 480), 'PH{c} highest chs [SNR]', tempcsnr, ch == self.cluster_size)
			minz, maxz = self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)]
			PlotCellsProfiles('PH{c}_H_cell_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], minz, max(maxz, 4800), 'PH{c} highest chs [ADC]', tempcadc, ch == self.cluster_size)

	def DoFinalAnalysis(self):
		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz
		self.GetCutsFromCutManager('all')
		self.GetCutsFromCutManager('good')
		self.GetVarzFromTranspGrid()
		self.DefineSatRegion(before=0, after=1)
		self.DoDeviceMaps('all', '', 'no_cuts')
		list_cuts = ['', 'no_neg', 'no_neg_no_sat']
		cells = 'good'
		for cut in list_cuts:
			self.DoDeviceMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut))
			self.DoCellMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut))

	def GetCutsFromCutManager(self, cells):
		self.noise_cuts[cells] = self.trans_grid.cuts_man.noise_cuts[cells]
		self.noise_friend_cuts[cells] = self.trans_grid.cuts_man.noise_friend_cuts[cells]
		self.in_transp_cluster = self.trans_grid.cuts_man.ConcatenateCuts(cut2=self.trans_grid.cuts_man.in_transp_cluster, cut1=self.trans_grid.cuts_man.transp_ev)
		self.noise_nc_cuts[cells] = self.trans_grid.cuts_man.noise_nc_cuts[cells]
		self.noise_nc_friend_cuts[cells] = self.trans_grid.cuts_man.noise_nc_friend_cuts[cells]

		self.ph_adc_ch_cuts[cells] = self.trans_grid.cuts_man.ph_adc_ch_cuts[cells]
		self.ph_snr_ch_cuts[cells] = self.trans_grid.cuts_man.ph_snr_ch_cuts[cells]
		self.ph_adc_h_cuts[cells] = self.trans_grid.cuts_man.ph_adc_h_cuts[cells]
		self.ph_snr_h_cuts[cells] = self.trans_grid.cuts_man.ph_snr_h_cuts[cells]
		self.phN_adc_ch_cuts[cells] = self.trans_grid.cuts_man.phN_adc_ch_cuts[cells]
		self.phN_snr_ch_cuts[cells] = self.trans_grid.cuts_man.phN_snr_ch_cuts[cells]
		self.phN_adc_h_cuts[cells] = self.trans_grid.cuts_man.phN_adc_h_cuts[cells]
		self.phN_snr_h_cuts[cells] = self.trans_grid.cuts_man.phN_snr_h_cuts[cells]

		self.sat_adc_N_ch_cut = self.trans_grid.cuts_man.sat_adc_N_ch
		self.sat_adc_N_h_cut = self.trans_grid.cuts_man.sat_adc_N_h
		self.not_sat_adc_N_ch_cut = self.trans_grid.cuts_man.not_sat_adc_N_ch
		self.not_sat_adc_N_h_cut = self.trans_grid.cuts_man.not_sat_adc_N_h

		self.sat_evts_region = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region = self.trans_grid.cuts_man.not_sat_evts_region

		self.neg_snr_phN_ch = self.trans_grid.cuts_man.neg_snr_phN_ch
		self.neg_snr_phN_h = self.trans_grid.cuts_man.neg_snr_phN_h
		self.not_neg_snr_phN_ch = self.trans_grid.cuts_man.not_neg_snr_phN_ch
		self.not_neg_snr_phN_h = self.trans_grid.cuts_man.not_neg_snr_phN_h

		self.neg_adc_phN_ch = self.trans_grid.cuts_man.neg_adc_phN_ch
		self.neg_adc_phN_h = self.trans_grid.cuts_man.neg_adc_phN_h
		self.not_neg_adc_phN_ch = self.trans_grid.cuts_man.not_neg_adc_phN_ch
		self.not_neg_adc_phN_h = self.trans_grid.cuts_man.not_neg_adc_phN_h

	def GetVarzFromTranspGrid(self):
		self.noise_varz = self.trans_grid.noise_varz
		self.ph_adc_ch_varz = self.trans_grid.ph_adc_ch_varz
		self.ph_snr_ch_varz = self.trans_grid.ph_snr_ch_varz
		self.ph_adc_h_varz = self.trans_grid.ph_adc_h_varz
		self.ph_snr_h_varz = self.trans_grid.ph_snr_h_varz
		self.phN_adc_ch_varz = self.trans_grid.phN_adc_ch_varz
		self.phN_snr_ch_varz = self.trans_grid.phN_snr_ch_varz
		self.phN_adc_h_varz = self.trans_grid.phN_adc_h_varz
		self.phN_snr_h_varz = self.trans_grid.phN_snr_h_varz

	def DefineSatRegion(self, before=0, after=1):
		if self.trans_grid.trans_tree.GetFriend('satRegions'):
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('satRegions'))
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=after, skipBefore=before)

if __name__ == '__main__':
	c = FinalAnalysis(None, 0, 0)