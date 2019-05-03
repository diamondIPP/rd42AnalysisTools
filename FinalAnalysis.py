#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from NoiseAnalysis import NoiseAnalysis
from CenterCellAnalysis import CenterCellAnalysis
from Utils import *


color_index = 10000

class FinalAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana=None, center_reg_ana=None):
		self.window_shift = 2
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
		self.center_reg_ana = center_reg_ana
		self.cell_dists = np.arange(0.1, 1, 0.025)

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz

		self.phN_chs_var = self.trans_grid.GetPHNChsVar
		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

		self.sat_evts_region_cut = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region_cut = self.trans_grid.cuts_man.not_sat_evts_region

		self.neg_ch_cut = self.trans_grid.cuts_man.GetNegPHChCut
		self.negN_chs_cut = self.trans_grid.cuts_man.GetNegPHNChsCut
		self.not_neg_ch_cut = self.trans_grid.cuts_man.GetNotNegPHChCut
		self.not_negN_chs_cut = self.trans_grid.cuts_man.GetNotNegPHNChsCut

		self.analysis_cummulative_ch = np.arange(1, self.num_strips + 1)

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def GetCut(self, cuts, typ, isFriend, ch=0, chtype='Ch'):
		chs = self.cluster_size if ch == 0 else ch
		return self.not_negN_chs_cut(chs, chtype, typ == 'snr', isFriend) if cuts == 'no_neg' else self.trans_grid.cuts_man.AndCuts([self.not_negN_chs_cut(chs, chtype, typ == 'snr', isFriend), self.not_sat_evts_region_cut]) if cuts == 'no_neg_no_sat' else '(1)'

	def DoDeviceMaps(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		def DrawProfile2D(name, varz, varzname, cut, getOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondMap(name, varz, varzname, cells, cut, 'prof colz')
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(self.trans_grid.row_cell_info_diamond['0'] - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_cell_info_diamond['up'] + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)
			if getOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.trans_grid.histo['hit_map_' + name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_' + name].GetYaxis().SetRangeUser(self.trans_grid.row_cell_info_diamond['0'] - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_cell_info_diamond['up'] + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
				self.trans_grid.DrawTCutGs('hit_map_' + name, 'diamond')
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_' + name)
				self.PosCanvas('hit_map_' + name)

		tempc = self.GetCut(cuts, typ, isFriend)
		for ch in self.analysis_cummulative_ch:
			hname = 'PH{c}_Ch_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			DrawProfile2D(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempc, ch == self.analysis_cummulative_ch[0])
			if ch != self.cluster_size:
				hname = 'PH{c}_H_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				DrawProfile2D(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempc, False)

	def DoStripHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_cell_info_diamond['height']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'
		def Draw2DHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			histo_limits = Get1DLimits(zmin, zmax, deltay)
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.PosCanvas(name)

		def Draw1DHistogram(name, cuts):
			self.trans_grid.DrawHisto1D(name, minx, maxx, deltax, xvar, xname, cuts)
			maxbin = self.trans_grid.histo[name].GetMaximumBin()
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(0, self.trans_grid.histo[name].GetBinContent(maxbin) + self.trans_grid.histo[name].GetBinError(maxbin))
			self.PosCanvas(name)

		tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.GetCut(cuts, typ, isFriend), cells=cells)

		for ch in self.analysis_cummulative_ch:
			minz, maxz = min(0, self.trans_grid.minz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]), self.trans_grid.maxz['PH{c}_Ch_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_Ch_Vs_strip_location_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), tempc, typ)
			minz, maxz = min(0, self.trans_grid.minz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]), self.trans_grid.maxz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
			if ch != self.cluster_size:
				hname = 'PH{c}_H_Vs_strip_location_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_{t}_{s}'.format(c=ch, t=typ.lower(), s=suffix, v=self.trans_grid.noise_friend_buffer)
				Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempc, typ)
		Draw1DHistogram('strip_location_{t}_{s}'.format(s=suffix, t=typ.lower()), tempc)

	def DoPH2DHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.GetCut(cuts, typ, isFriend), cells=cells)
		for ch in self.analysis_cummulative_ch:
			nameh = 'PH{c}_H_Vs_hit_channel_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_hit_channel_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			xmin, xmax, deltax = self.trans_grid.ch_ini, self.trans_grid.ch_end, 1
			ymin, ymax, deltay = (self.trans_grid.phmin, self.trans_grid.phmax, float(self.trans_grid.phmax) / self.trans_grid.phbins) if typ == 'adc' else (self.trans_grid.phmin, self.trans_grid.phmax / 10., float(self.trans_grid.phmax) / self.trans_grid.phbins / 10.)
			vary = self.phN_chs_var(ch, 'H', typ == 'snr', isFriend)
			self.trans_grid.DrawHisto2D(nameh, xmin, xmax, deltax, 'dia pred hit ch', ymin, ymax, deltay, 'PH{c} highest chs'.format(c=ch), 'diaChannels[int(TMath::Floor(diaChXPred+0.5))]', vary, tempc)
			self.PosCanvas(nameh)
			nameh = 'PH{c}_H_Vs_event_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_event_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			xmin, xmax, deltax = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event'), 100 * self.delta_ev
			self.trans_grid.DrawHisto2D(nameh, xmin, xmax, deltax, 'event', ymin, ymax, deltay, 'PH{c} highest chs'.format(c=ch), 'event', vary, tempc)
			self.PosCanvas(nameh)

			namehp = 'PH{c}_H_mean_Vs_event_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_mean_Vs_event_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			self.trans_grid.histo[namehp] = self.trans_grid.histo[nameh].ProfileX('h_' + namehp)
			self.trans_grid.histo[namehp].SetTitle('h_' + namehp)
			self.trans_grid.canvas[namehp] = ro.TCanvas('c_' + namehp, 'c_' + namehp, 1)
			self.trans_grid.histo[namehp].Draw('e hist')
			SetDefault1DStats(self.trans_grid.histo[namehp], y1=0.15, y2=0.45)
			SetDefault1DCanvasSettings(self.trans_grid.canvas[namehp])
			ro.gPad.Update()
			self.trans_grid.FitPol(namehp, 1)
			self.PosCanvas(namehp)

	def DoEfficiencyPlots(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		def DrawEfficiencyGraphs(name, var, cells, cuts, typ='adc', show_only_95=True):
			xmin, xmax, deltax = (self.trans_grid.phmin, self.trans_grid.phmax, 50) if typ == 'adc' else (self.trans_grid.phmin / 10., self.trans_grid.phmax / 10., 5)
			ymin = 0.95 if show_only_95 else 0
			self.trans_grid.DrawEfficiencyGraph(name, var, cells, cuts, xmin, xmax, deltax, typ, ymin)
			self.PosCanvas(name)

		def DrawEfficiencyCellMaps(name, var, cells, cuts, ncuts=[1, 2, 5, 10, 20], typ='adc'):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name + '_h0_', var=var, cells=cells, cuts=cuts, plot_option='prof colz goff')
			self.trans_grid.GetOccupancyFromProfile(name + '_h0_', 'colz goff')
			typ_cut = 'sigma' if typ == 'snr' else 'adc'
			for ns in ncuts:
				tempcut = '({v}>={th})'.format(v=var, th=ns)
				cut = tempcut if cuts == '' else self.trans_grid.cuts_man.AndCuts([cuts, tempcut])
				# cut = tempcut if cuts == '' else self.trans_grid.cuts_man.ConcatenateCuts(cuts, tempcut)
				self.trans_grid.DrawProfile2DDiamondCellOverlay(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), var=var, cells=cells, cuts=cut, plot_option='prof colz goff')
				self.trans_grid.GetOccupancyFromProfile(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 'colz goff')
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)] = self.trans_grid.histo['hit_map_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Clone('h_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].SetTitle('h_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Sumw2()
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Divide(self.trans_grid.histo['hit_map_' + name + '_h0_'])
				self.trans_grid.canvas[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)] = ro.TCanvas('c_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 'c_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 1)
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].GetZaxis().SetTitle('Efficiency')
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].GetZaxis().SetRangeUser(0.9, 1)
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Draw('colz')
				SetDefault2DStats(self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)])
				self.PosCanvas(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))

		tempc = self.GetCut(cuts, typ, isFriend)
		for ch in self.analysis_cummulative_ch:
			hname = 'PH{c}_H_Efficiency_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_Efficiency_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			ncuts = [10, 50, 100, 200, 500] if typ == 'adc' else [1, 5, 10, 20, 50]
			DrawEfficiencyCellMaps(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), cells, tempc, ncuts=ncuts, typ=typ)
			hname = 'Eff_PH{c}_H_Vs_Threshold_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else'Eff_PH{c}_H_buffer_{v}_Vs_Threshold_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			DrawEfficiencyGraphs(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), cells, tempc, typ)

	def DoPHHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		def DrawHisto(name, xmin, xmax, deltax, varz, varname, cuts):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varname, xmin=xmin, xmax=xmax, deltax=deltax)
			self.PosCanvas(name)

		tempc = self.GetCut(cuts, typ, isFriend)
		xmin, xmax = (self.trans_grid.phmin, self.trans_grid.phmax) if typ == 'adc' else (self.trans_grid.phmin / 10., self.trans_grid.phmax / 10.)
		deltax = float(xmax - xmin) / self.trans_grid.phbins
		for ch in self.analysis_cummulative_ch:
			for chtype in ['Ch', 'H']:
				if ch != self.cluster_size or chtype != 'H':
					hname = 'PH{c}_{ct}_{t}_{s}'.format(c=ch, ct=chtype, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_{ct}_buffer_{v}_{t}_{s}'.format(c=ch, ct=chtype, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
					DrawHisto(hname, xmin, xmax, deltax, self.phN_chs_var(ch, chtype, typ == 'snr', isFriend), 'PH{c} {ct} chs [{t}]'.format(c=ch, ct='cluster' if chtype == 'Ch' else 'highest', t=typ.upper()), tempc)
					self.trans_grid.FitLanGaus(hname)

	def DoCenterPHStudies(self, cells, cuts='', suffix='no_cuts', typ='adc', do_sat=True, do_all_plots=False, isFriend=False):
		self.center_reg_ana.w, self.center_reg_ana.window_shift = self.w, self.window_shift
		tempc = self.GetCut(cuts, typ, isFriend)
		self.center_reg_ana.DoCenterRegionStudies(self.cell_dists, cells, do_all_plots, suffix, tempc, typ, do_sat, isFriend)
		self.w, self.window_shift = self.center_reg_ana.w, self.center_reg_ana.window_shift

	def DoCellMaps(self, cells, cuts='', suffix='no_cuts', typ='adc', isFriend=False):
		def PlotCellsProfiles(name, varz, zmin, zmax, varname, cut, doOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname)
			self.trans_grid.profile[name].SetMinimum(min(0, zmin))
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.PosCanvas(name)
			if doOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.PosCanvas('hit_map_' + name)

		tempc = self.GetCut(cuts, typ, isFriend)
		for ch in self.analysis_cummulative_ch:
			for chtype in ['Ch', 'H']:
				if ch != self.cluster_size or chtype != 'H':
					minz, maxz = self.minz['PH{c}_{ct}_{t}'.format(c=ch, t=typ.lower(), ct=chtype)], max(self.minz['PH{c}_{ct}_{t}'.format(c=ch, t=typ.lower(), ct=chtype)], self.trans_grid.phmax if typ == 'adc' else self.trans_grid.phmax / 10.)
					hname = 'PH{c}_{ct}_cell_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), ct=chtype) if not isFriend else 'PH{c}_{ct}_buffer_{v}_cell_map_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer, ct=chtype)
					PlotCellsProfiles(hname, self.phN_chs_var(ch, chtype, typ == 'snr', isFriend), minz, maxz, 'PH{c} {ct} chs [{t}]'.format(c=ch, t=typ.upper(), ct='cluster' if chtype == 'Ch' else 'highest'), tempc)

	def DoFinalAnalysis(self, typ='adc', cummulative_chs=None, isFriend=False):
		if cummulative_chs:
			self.analysis_cummulative_ch = cummulative_chs
		self.DefineSatRegion(before=0, after=1)
		self.DoDeviceMaps('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoCellMaps('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoStripHistograms('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoPH2DHistograms('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoEfficiencyPlots('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoPHHistograms('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		self.DoCenterPHStudies('all', '', 'no_cuts', typ=typ, isFriend=isFriend)
		list_cuts = ['', 'no_neg', 'no_neg_no_sat']
		cells = 'good'
		for cut in list_cuts:
			self.DoDeviceMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoCellMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoStripHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoPH2DHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoEfficiencyPlots(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoPHHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, isFriend=isFriend)
			self.DoCenterPHStudies(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, do_sat=(cut != 'no_neg_no_sat'), isFriend=isFriend)

	def DefineSatRegion(self, before=0, after=1):
		if self.trans_grid.trans_tree.GetFriend('satRegions'):
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('satRegions'))
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=after, skipBefore=before)

if __name__ == '__main__':
	c = FinalAnalysis(None, 0, 0)
