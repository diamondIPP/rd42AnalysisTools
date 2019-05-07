#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class SaturationAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana=None):
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

		self.suffix = GetSuffixDictionary(self.trans_grid)

		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz

		self.phN_chs_var = self.trans_grid.GetPHNChsVar
		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.sat_N_ch_cut = self.trans_grid.cuts_man.GetSatNChsCut

		self.sat_evts_region_cut = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region_cut = self.trans_grid.cuts_man.not_sat_evts_region

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

		self.strips_for_analysis = np.arange(self.cluster_size, dtype='i4')

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def DoProfileMaps(self, before=0, after=1, typ='adc', isFriend=False):
		"""
		Makes the 2D profile maps of events whose cluster have a saturated channel. The x and y axis are the predicted hit position, and the z axis is the PH of the N highest signals, where N is defined by the num_strips variable.
		:param before: number of events to take into account before the saturation event without including the saturated event
		:param after: number of events to take into account after the saturation event including the saturated event
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		self.DefineRegion(before=before, after=after)

		lims_dic = self.trans_grid.GetDiamondMapPlotLimits()
		xmin, xmax, ymin, ymax, xname, yname = lims_dic['xmin'], lims_dic['xmax'], lims_dic['ymin'], lims_dic['ymax'], 'pred dia hit ch', 'sil pred dia hit in Y [#mum]'
		deltax, deltay = float(xmax - xmin) / lims_dic['binsx'], float(ymax - ymin) / lims_dic['binsy']
		yup_max = max(self.trans_grid.row_cell_info_diamond['up_odd'], self.trans_grid.row_cell_info_diamond['up_even'])
		ylow_min = min(self.trans_grid.row_cell_info_diamond['0_even'], self.trans_grid.row_cell_info_diamond['0_odd'])

		def DrawProfile2D(name, varz, zmin, zmax, varzname, cut, xdelt=deltax, namex=xname, varx='diaChXPred', getOccupancy=False):
			self.trans_grid.DrawProfile2D(name, xmin, xmax, xdelt, namex, ymin, ymax, deltay, yname, varx, 'diaChYPred', varz, varzname, cut)
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(ylow_min - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), yup_max + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.profile[name].SetMinimum(zmin)
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)
			if getOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.trans_grid.histo['hit_map_' + name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_' + name].GetYaxis().SetRangeUser(ylow_min - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), yup_max + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
				self.trans_grid.DrawTCutGs('hit_map_' + name, 'diamond')
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_' + name)
				self.PosCanvas('hit_map_' + name)


		tempc = self.sat_N_ch_cut(self.num_strips, 'H')
		minz, maxz = min(self.minz['PH{n}_H_{t}'.format(n=self.num_strips, t=typ.lower())], 0), self.maxz['PH{n}_H_{t}'.format(n=self.num_strips, t=typ.lower())]
		hname = 'PH{n}_H_pred_hit_sat_events_map_{b}_before_{a}_after_{t}'.format(n=self.num_strips, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{n}_H_buffer_{v}_pred_hit_sat_events_map_{b}_before_{a}_after_{t}'.format(n=self.num_strips, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
		DrawProfile2D(hname, self.phN_chs_var(self.num_strips, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempc, getOccupancy=True)
		tempc = self.sat_N_ch_cut(1, 'H')
		minz, maxz = min(self.minz['PH1_H_{t}'.format(t=typ.lower())], 0), self.maxz['PH1_H_{t}'.format(t=typ.lower())]
		hname = 'PH_sat_ch_map_{b}_before_{a}_after_{t}'.format(b=before, a=after, t=typ.lower()) if not isFriend else 'PH_buffer_{v}_sat_ch_map_{b}_before_{a}_after_{t}'.format(b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
		DrawProfile2D(hname, self.phN_chs_var(1, 'H', typ == 'snr', isFriend), minz, maxz, 'PH sat ch [{t}]'.format(t=typ.upper()), tempc, 1, 'highest ch', 'clusterChannelHighest1')

	def Do1DPHHistos(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		"""
		Makes the PH distribution of several channels ordered either by proximity from the hit position or by highest signal, for saturated clusters
		:param cells: Which cells to take into account for the profile
		:param before: number of events to take into account before the saturation event without including the saturated event
		:param after: number of events to take into account after the saturation event including the saturated event
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		self.DefineRegion(before=before, after=after)
		suffix = self.suffix[cells]
		def DrawPHHisto(name, varz, varzname, cuts):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varzname, typ=typ)
			self.PosCanvas(name)
		tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		for ch in (self.strips_for_analysis + 1):
			hname = 'PH{c}_Ch_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			DrawPHHisto(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat)
			hname = 'PH{c}_Ch_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			DrawPHHisto(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat)
			if ch != self.cluster_size:
				hname = 'PH{c}_H_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				DrawPHHisto(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat)
				hname = 'PH{c}_H_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				DrawPHHisto(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat)

	def DoCellMaps(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		"""
		Creates cell profile maps of overlaid cells for saturated clusters. The Z axis is the cumulative PH of different channels
		:param cells: Which cells to take into account for the profile
		:param before: number of events to take into account before the saturation event without including the saturated event
		:param after: number of events to take into account after the saturation event including the saturated event
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		self.DefineRegion(before=before, after=after)
		def PlotCellsProfiles(name, varz, zmin, zmax, varname, cut, doOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname)
			self.trans_grid.profile[name].SetMinimum(zmin)
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.PosCanvas(name)
			if doOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.PosCanvas('hit_map_' + name)

		suffix = self.suffix[cells]
		for ch in (self.strips_for_analysis + 1):
			tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=ch), isFriend)])
			tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut,  self.ph_cuts('PH{c}_Ch'.format(c=ch), isFriend)])
			minz, maxz = min(self.minz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())], 0), self.maxz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_Ch_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			PlotCellsProfiles(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat, ch == 1)
			hname = 'PH{c}_Ch_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			PlotCellsProfiles(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat, ch == 1)
			tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_H'.format(c=ch), isFriend)])
			tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, self.ph_cuts('PH{c}_H'.format(c=ch), isFriend)])
			minz, maxz = min(self.minz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())], 0), self.maxz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
			if ch != self.cluster_size:
				hname = 'PH{c}_H_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				PlotCellsProfiles(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat)
				hname = 'PH{c}_H_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				PlotCellsProfiles(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat)

	def PlotStripHistograms(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		"""
		Creates cell histograms indicating the position where the saturated cluster hit the hit strip (ch0) and the PH for different channel configurations
		:param cells: Which cells to take into account for the profile
		:param before: number of events to take into account before the saturation event without including the saturated event
		:param after: number of events to take into account after the saturation event including the saturated event
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		self.DefineRegion(before=before, after=after)
		minx, maxx, deltax, xname, xvar = -self.trans_grid.row_cell_info_diamond['width'] / (2.0 *self.trans_grid.col_pitch), self.trans_grid.row_cell_info_diamond['width'] / (2.0 *self.trans_grid.col_pitch), self.trans_grid.cell_resolution / float(self.trans_grid.row_cell_info_diamond['width']), 'dia pred. strip hit pos', 'x0'

		def Draw2DHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			histo_limits = Get1DLimits(zmin, zmax, 4 * deltay)
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, tempc)
			self.PosCanvas(name)

		def Draw1DHistogram(name, cuts):
			tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
			self.trans_grid.DrawHisto1D(name, minx, maxx, deltax, xvar, xname, tempc)
			maxbin = self.trans_grid.histo[name].GetMaximumBin()
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(0, self.trans_grid.histo[name].GetBinContent(maxbin) + self.trans_grid.histo[name].GetBinError(maxbin))
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''

		tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		for ch in (self.strips_for_analysis + 1):
			minz, maxz = self.trans_grid.minz['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())], self.trans_grid.maxz['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())]
			hname = 'PH_Ch{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch-1, t=typ.upper()), self.ph_ch_var(ch - 1, 'Ch', typ == 'snr', isFriend), tempcsat)
			hname = 'PH_Ch{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch-1, t=typ.upper()), self.ph_ch_var(ch - 1, 'Ch', typ == 'snr', isFriend), tempcnosat)
			minz, maxz = self.trans_grid.minz['PH_H{c}_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz['PH_H{c}_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH_H{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH_H{c}_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(c=ch, t=typ.upper(), sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_ch_var(ch, 'H', typ == 'snr', isFriend), tempcsat)
			hname = 'PH_H{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH_H{c}_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(c=ch, t=typ.upper(), sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_ch_var(ch, 'H', typ == 'snr', isFriend), tempcnosat)
			minz, maxz = self.trans_grid.minz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_Ch_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), tempcsat)
			hname = 'PH{c}_Ch_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), tempcnosat)
			if ch != self.cluster_size or ch != 1:
				minz, maxz = self.trans_grid.minz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
				hname = 'PH{c}_H_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempcsat)
				hname = 'PH{c}_H_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempcnosat)

		Draw1DHistogram('strip_location_not_sat_events_{b}_before_{a}_after_{s}'.format(b=before, a=after, s=suffix), tempcnosat)
		Draw1DHistogram('strip_location_sat_events_{b}_before_{a}_after_{s}'.format(b=before, a=after, s=suffix), tempcsat)

	def DoSaturationAnalysis(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		self.SetStripsForAnalysis()
		self.DoProfileMaps(before=0, after=1, typ=typ, isFriend=isFriend)
		self.DoCellMaps(cells, before=0, after=1, typ=typ, isFriend=isFriend)
		self.PlotStripHistograms(cells, before=0, after=1, typ=typ, isFriend=isFriend)
		self.Do1DPHHistos(cells, before=0, after=1, typ=typ, isFriend=isFriend)
		if before != 0 or after != 1:
			self.PlotStripHistograms(cells, before=before, after=after, typ=typ, isFriend=isFriend)
			self.Do1DPHHistos(cells, before=before, after=after, typ=typ, isFriend=isFriend)

	def DefineRegion(self, before=0, after=1):
		if self.trans_grid.trans_tree.GetFriend('satRegions'):
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('satRegions'))
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=after, skipBefore=before)

	def SetStripsForAnalysis(self, arr=None):
		if arr:
			self.strips_for_analysis = arr
		else:
			self.strips_for_analysis = np.arange(self.num_strips, dtype='i4')
			self.strips_for_analysis = self.strips_for_analysis if self.cluster_size > 3 else np.unique(np.append(self.strips_for_analysis, self.cluster_size - 1))


if __name__ == '__main__':
	c = SaturationAnalysis(None, 0, 0)
