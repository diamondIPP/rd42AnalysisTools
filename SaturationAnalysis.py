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

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz

		self.phN_chs_var = self.trans_grid.GetPHNChsVar
		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.sat_N_ch_cut = self.trans_grid.cuts_man.GetSatNChsCut

		self.sat_evts_region_cut = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region_cut = self.trans_grid.cuts_man.not_sat_evts_region

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def DoProfileMaps(self, before=0, after=1, typ='adc', isFriend=False):
		self.DefineRegion(before=before, after=after)

		xmin, xmax, deltax, xname = self.trans_grid.ch_ini - 1.5, self.trans_grid.ch_end + 1.5, 1.0/self.trans_grid.bins_per_ch_x, 'pred dia hit ch',
		ymin, ymax, deltay, yname = self.trans_grid.row_info_diamond['0'] - RoundInt(float(self.trans_grid.row_info_diamond['0']) / self.trans_grid.row_info_diamond['pitch'], 'f8') * self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['0'] + (256 - RoundInt(float(self.trans_grid.row_info_diamond['0']) / self.trans_grid.row_info_diamond['pitch'], 'f8')) * self.trans_grid.row_info_diamond['pitch'], float(self.trans_grid.row_info_diamond['pitch'])/self.trans_grid.bins_per_ch_y, 'sil pred dia hit in Y [#mum]'

		def DrawProfile2D(name, varz, zmin, zmax, varzname, cut, xdelt=deltax, namex=xname, varx='diaChXPred', getOccupancy=False):
			self.trans_grid.DrawProfile2D(name, xmin, xmax, xdelt, namex, ymin, ymax, deltay, yname, varx, 'diaChYPred', varz, varzname, cut)
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.profile[name].SetMinimum(zmin)
			self.trans_grid.profile[name].SetMaximum(zmax)
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


		tempc = self.sat_N_ch_cut(self.num_strips, 'H')
		minz, maxz = min(self.minz['all']['PH{n}_H_{t}'.format(n=self.num_strips, t=typ.lower())], 0), self.maxz['all']['PH{n}_H_{t}'.format(n=self.num_strips, t=typ.lower())]
		hname = 'PH{n}_H_pred_hit_sat_events_map_{b}_before_{a}_after_{t}'.format(n=self.num_strips, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{n}_H_buffer_{v}_pred_hit_sat_events_map_{b}_before_{a}_after_{t}'.format(n=self.num_strips, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
		DrawProfile2D(hname, self.phN_chs_var(self.num_strips, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempc, getOccupancy=True)
		tempc = self.sat_N_ch_cut(1, 'H')
		minz, maxz = min(self.minz['all']['PH1_H_{t}'.format(t=typ.lower())], 0), self.maxz['all']['PH1_H_{t}'.format(t=typ.lower())]
		hname = 'PH_sat_ch_map_{b}_before_{a}_after_{t}'.format(b=before, a=after, t=typ.lower()) if not isFriend else 'PH_buffer_{v}_sat_ch_map_{b}_before_{a}_after_{t}'.format(b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
		DrawProfile2D(hname, self.phN_chs_var(1, 'H', typ == 'snr', isFriend), minz, maxz, 'PH sat ch [{t}]'.format(t=typ.upper()), tempc, 1, 'highest ch', 'clusterChannelHighest1')

	def Do1DPHHistos(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		self.DefineRegion(before=before, after=after)
		suffix = self.suffix[cells]
		def DrawPHHisto(name, varz, varzname, cuts):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varzname, typ=typ)
			self.PosCanvas(name)
		tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=self.cluster_size), isFriend)])
		for ch in xrange(1, self.cluster_size + 1):
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
		for ch in xrange(1, self.cluster_size + 1):
			tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_Ch'.format(c=ch), isFriend)])
			tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut,  self.ph_cuts('PH{c}_Ch'.format(c=ch), isFriend)])
			minz, maxz = min(self.minz[cells]['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())], 0), self.maxz[cells]['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_Ch_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			PlotCellsProfiles(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat, ch == 1)
			hname = 'PH{c}_Ch_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			PlotCellsProfiles(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat, ch == 1)
			tempcsat = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, self.ph_cuts('PH{c}_H'.format(c=ch), isFriend)])
			tempcnosat = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, self.ph_cuts('PH{c}_H'.format(c=ch), isFriend)])
			minz, maxz = min(self.minz[cells]['PH{c}_H_{t}'.format(c=ch, t=typ.lower())], 0), self.maxz[cells]['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
			if ch != self.cluster_size:
				hname = 'PH{c}_H_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_cell_map_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				PlotCellsProfiles(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcsat)
				hname = 'PH{c}_H_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_cell_map_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				PlotCellsProfiles(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempcnosat)

	def PlotStripHistograms(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
		self.DefineRegion(before=before, after=after)
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_info_diamond['pitch']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'

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
		for ch in xrange(1, self.cluster_size + 1):
			minz, maxz = self.trans_grid.minz[cells]['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())], self.trans_grid.maxz[cells]['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())]
			hname = 'PH_Ch{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch-1, t=typ.upper()), self.ph_ch_var(ch - 1, 'Ch', typ == 'snr', isFriend), tempcsat)
			hname = 'PH_Ch{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch-1, s=suffix, b=before, a=after, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch-1, t=typ.upper()), self.ph_ch_var(ch - 1, 'Ch', typ == 'snr', isFriend), tempcnosat)
			minz, maxz = self.trans_grid.minz[cells]['PH_H{c}_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz[cells]['PH_H{c}_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH_H{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH_H{c}_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(c=ch, t=typ.upper(), sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_ch_var(ch, 'H', typ == 'snr', isFriend), tempcsat)
			hname = 'PH_H{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH_H{c}_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(c=ch, t=typ.upper(), sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_ch_var(ch, 'H', typ == 'snr', isFriend), tempcnosat)
			minz, maxz = self.trans_grid.minz[cells]['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz[cells]['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_Ch_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), tempcsat)
			hname = 'PH{c}_Ch_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), tempcnosat)
			if ch != self.cluster_size or ch != 1:
				minz, maxz = self.trans_grid.minz[cells]['PH{c}_H_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz[cells]['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
				hname = 'PH{c}_H_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempcsat)
				hname = 'PH{c}_H_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, b=before, a=after, v=self.trans_grid.noise_friend_buffer)
				Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempcnosat)

		Draw1DHistogram('strip_location_not_sat_events_{b}_before_{a}_after_{s}'.format(b=before, a=after, s=suffix), tempcnosat)
		Draw1DHistogram('strip_location_sat_events_{b}_before_{a}_after_{s}'.format(b=before, a=after, s=suffix), tempcsat)

	def DoSaturationAnalysis(self, cells='all', before=0, after=1, typ='adc', isFriend=False):
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


if __name__ == '__main__':
	c = SaturationAnalysis(None, 0, 0)
