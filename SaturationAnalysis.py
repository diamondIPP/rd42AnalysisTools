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

	def DoProfileMaps(self, before=0, after=1):
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

		tempc = self.sat_adc_N_h_cut['2_H']
		minz, maxz = min(self.minz['all']['PH2_H_snr'], 0), self.maxz['all']['PH2_H_snr']
		DrawProfile2D('PH2_H_pred_hit_sat_events_map_{b}_before_{a}_after_snr'.format(b=before, a=after), self.phN_snr_h_varz['PH2_H'], minz, maxz, 'PH2 highest chs [SNR]', tempc, getOccupancy=True)
		minz, maxz = min(self.minz['all']['PH2_H_adc'], 0), self.maxz['all']['PH2_H_adc']
		DrawProfile2D('PH2_H_pred_hit_sat_events_map_{b}_before_{a}_after_adc'.format(b=before, a=after), self.phN_adc_h_varz['PH2_H'], minz, maxz, 'PH2 highest chs [ADC]', tempc, getOccupancy=True)
		tempc = self.sat_adc_N_h_cut['1_H']
		minz, maxz = min(self.minz['all']['PH1_H_snr'], 0), self.maxz['all']['PH1_H_snr']
		DrawProfile2D('PH_sat_ch_map_{b}_before_{a}_after_snr'.format(b=before, a=after), self.phN_snr_h_varz['PH1_H'], minz, maxz, 'PH sat ch [snr]', tempc, 1, 'highest ch', 'clusterChannelHighest1')
		minz, maxz = min(self.minz['all']['PH1_H_adc'], 0), self.maxz['all']['PH1_H_adc']
		DrawProfile2D('PH_sat_ch_map_{b}_before_{a}_after_adc'.format(b=before, a=after), self.phN_adc_h_varz['PH1_H'], minz, maxz, 'PH sat ch [ADC]', tempc, 1, 'highest ch', 'clusterChannelHighest1')

	def Do1DPHHistos(self, cells='all', before=0, after=1):
		self.DefineRegion(before=before, after=after)
		suffix = self.suffix[cells]
		noise_name0 = 'signal_noise_{s}_{t}'.format(s=suffix, t='adc')
		if not noise_name0 in self.trans_grid.histo.keys():
			self.noise_ana.PlotNoiseNotInCluster(cells)
		sigm = self.trans_grid.histo[noise_name0].GetRMS()
		def DrawPHHisto(name, varz, varzname, cuts, xmin=10000, xmax=-10000, deltax=-1):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varzname, xmin=xmin, xmax=xmax, deltax=deltax)
			self.PosCanvas(name)
		tempcsat = self.trans_grid.cuts_man.ConcatenateCuts(self.sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=self.cluster_size)])
		tempcnosat = self.trans_grid.cuts_man.ConcatenateCuts(self.not_sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=self.cluster_size)])
		minsnr, maxsnr = int(RoundInt(self.trans_grid.phmin / sigm)), int(RoundInt(self.trans_grid.phmax / sigm))
		deltsnr = float(maxsnr - minsnr) / float(self.trans_grid.phbins)
		for ch in xrange(1, self.cluster_size + 1):
			tempc = tempcsat
			DrawPHHisto('PH{c}_Ch_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} cluster chs [SNR]', tempc, minsnr, maxsnr, deltsnr)
			tempc = tempcnosat
			DrawPHHisto('PH{c}_Ch_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} cluster chs [ADC]', tempc)

			tempc = tempcsat
			DrawPHHisto('PH{c}_H_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [SNR]', tempc, minsnr, maxsnr, deltsnr)
			tempc = tempcnosat
			DrawPHHisto('PH{c}_H_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [ADC]', tempc)

	def DoCellMaps(self, cells='all', before=0, after=1):
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
			tempcsat = self.trans_grid.cuts_man.ConcatenateCuts(self.sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=ch)])
			tempcnosat = self.trans_grid.cuts_man.ConcatenateCuts(self.not_sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=ch)])
			minzadc, maxzadc = min(self.minz[cells]['PH{c}_Ch_adc'.format(c=ch)], 0), self.maxz[cells]['PH{c}_Ch_adc'.format(c=ch)]
			minzsnr, maxzsnr = min(self.minz[cells]['PH{c}_Ch_snr'.format(c=ch)], 0), self.maxz[cells]['PH{c}_Ch_snr'.format(c=ch)]
			PlotCellsProfiles('PH{c}_Ch_cell_map_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], minzsnr, maxzsnr, 'PH{c} cluster chs [SNR]'.format(c=ch), tempcsat)
			PlotCellsProfiles('PH{c}_Ch_cell_map_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], minzadc, maxzadc, 'PH{c} cluster chs [ADC]'.format(c=ch), tempcsat, True)
			PlotCellsProfiles('PH{c}_Ch_cell_map_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], minzsnr, maxzsnr, 'PH{c} cluster chs [SNR]'.format(c=ch), tempcnosat)
			PlotCellsProfiles('PH{c}_Ch_cell_map_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], minzadc, maxzadc, 'PH{c} cluster chs [ADC]'.format(c=ch), tempcnosat, True)
			tempcsat = self.trans_grid.cuts_man.ConcatenateCuts(self.sat_evts_region, self.phN_snr_h_cuts[cells]['PH{c}_H'.format(c=ch)])
			tempcnosat = self.trans_grid.cuts_man.ConcatenateCuts(self.not_sat_evts_region, self.phN_snr_h_cuts[cells]['PH{c}_H'.format(c=ch)])
			minzadc, maxzadc = min(self.minz[cells]['PH{c}_H_adc'.format(c=ch)], 0), self.maxz[cells]['PH{c}_H_adc'.format(c=ch)]
			minzsnr, maxzsnr = min(self.minz[cells]['PH{c}_H_snr'.format(c=ch)], 0), self.maxz[cells]['PH{c}_H_snr'.format(c=ch)]
			PlotCellsProfiles('PH{c}_H_cell_map_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], minzsnr, maxzsnr, 'PH{c} highest chs [SNR]'.format(c=ch), tempcsat)
			PlotCellsProfiles('PH{c}_H_cell_map_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], minzadc, maxzadc, 'PH{c} highest chs [ADC]'.format(c=ch), tempcsat, True)
			PlotCellsProfiles('PH{c}_H_cell_map_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], minzsnr, maxzsnr, 'PH{c} highest chs [SNR]'.format(c=ch), tempcnosat)
			PlotCellsProfiles('PH{c}_H_cell_map_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], minzadc, maxzadc, 'PH{c} highest chs [ADC]'.format(c=ch), tempcnosat, True)

	def PlotStripHistograms(self, cells='all', before=0, after=1):
		self.DefineRegion(before=before, after=after)
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_info_diamond['pitch']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'

		def Draw2DHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			histo_limits = Get1DLimits(zmin, zmax, 4 * self.delta_adc) if typ == 'adc' else Get1DLimits(zmin, zmax, 4 * self.delta_snr)
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.PosCanvas(name)

		def Draw1DHistogram(name, cuts):
			self.trans_grid.DrawHisto1D(name, minx, maxx, deltax, xvar, xname, cuts)
			maxbin = self.trans_grid.histo[name].GetMaximumBin()
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(0, self.trans_grid.histo[name].GetBinContent(maxbin) + self.trans_grid.histo[name].GetBinError(maxbin))
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''

		tempcsat = self.trans_grid.cuts_man.ConcatenateCuts(self.sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=self.cluster_size)])
		tempcnosat = self.trans_grid.cuts_man.ConcatenateCuts(self.not_sat_evts_region, self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=self.cluster_size)])
		for ch in xrange(1, self.cluster_size + 1):
			minzsnr, maxzsnr = self.trans_grid.minz[cells]['PH_Ch{c}_snr'.format(c=ch-1)], self.trans_grid.maxz[cells]['PH_Ch{c}_snr'.format(c=ch-1)]
			minzadc, maxzadc = self.trans_grid.minz[cells]['PH_Ch{c}_adc'.format(c=ch-1)], self.trans_grid.maxz[cells]['PH_Ch{c}_adc'.format(c=ch-1)]
			Draw2DHistogram('PH_Ch{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch-1, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH cluster ch{c} [SNR]'.format(c=ch-1), self.ph_snr_ch_varz['PH_Ch{c}'.format(c=ch-1)], tempcsat)
			Draw2DHistogram('PH_Ch{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch-1, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH cluster ch{c} [ADC]'.format(c=ch-1), self.ph_adc_ch_varz['PH_Ch{c}'.format(c=ch-1)], tempcsat)
			Draw2DHistogram('PH_Ch{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch-1, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH cluster ch{c} [SNR]'.format(c=ch-1), self.ph_snr_ch_varz['PH_Ch{c}'.format(c=ch-1)], tempcnosat)
			Draw2DHistogram('PH_Ch{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch-1, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH cluster ch{c} [ADC]'.format(c=ch-1), self.ph_adc_ch_varz['PH_Ch{c}'.format(c=ch-1)], tempcnosat)

			minzsnr, maxzsnr = self.trans_grid.minz[cells]['PH_H{c}_snr'.format(c=ch)], self.trans_grid.maxz[cells]['PH_H{c}_snr'.format(c=ch)]
			minzadc, maxzadc = self.trans_grid.minz[cells]['PH_H{c}_adc'.format(c=ch)], self.trans_grid.maxz[cells]['PH_H{c}_adc'.format(c=ch)]
			Draw2DHistogram('PH_H{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH highest {c}{sf} ch [SNR]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_snr_h_varz['PH_H{c}'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH_H{c}_Vs_strip_location_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH highest {c}{sf} ch [ADC]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_adc_h_varz['PH_H{c}'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH_H{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH highest {c}{sf} ch [SNR]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_snr_h_varz['PH_H{c}'.format(c=ch)], tempcnosat)
			Draw2DHistogram('PH_H{c}_Vs_strip_location_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH highest {c}{sf} ch [ADC]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_adc_h_varz['PH_H{c}'.format(c=ch)], tempcnosat)

			minzsnr, maxzsnr = self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)], self.trans_grid.maxz[cells]['PH{c}_Ch_snr'.format(c=ch)]
			minzadc, maxzadc = self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)], self.trans_grid.maxz[cells]['PH{c}_Ch_adc'.format(c=ch)]
			Draw2DHistogram('PH{c}_Ch_Vs_strip_location_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH{c} cluster chs [SNR]'.format(c=ch), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH{c}_Ch_Vs_strip_location_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH{c} cluster chs [ADC]'.format(c=ch), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH{c}_Ch_Vs_strip_location_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH{c} cluster chs [SNR]'.format(c=ch), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], tempcnosat)
			Draw2DHistogram('PH{c}_Ch_Vs_strip_location_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH{c} cluster chs [ADC]'.format(c=ch), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], tempcnosat)

			minzsnr, maxzsnr = self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)], self.trans_grid.maxz[cells]['PH{c}_H_snr'.format(c=ch)]
			minzadc, maxzadc = self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)], self.trans_grid.maxz[cells]['PH{c}_H_adc'.format(c=ch)]
			Draw2DHistogram('PH{c}_H_Vs_strip_location_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH{c} highest chs [SNR]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH{c}_H_Vs_strip_location_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH{c} highest chs [ADC]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], tempcsat)
			Draw2DHistogram('PH{c}_H_Vs_strip_location_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(c=ch, s=suffix, b=before, a=after), minzsnr, maxzsnr, 'PH{c} highest chs [SNR]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], tempcnosat)
			Draw2DHistogram('PH{c}_H_Vs_strip_location_not_sat_events_{b}_before_{a}_after_adc_{s}'.format(c=ch, s=suffix, b=before, a=after), minzadc, maxzadc, 'PH{c} highest chs [ADC]'.format(c=ch, sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], tempcnosat)

		Draw1DHistogram('strip_location_sat_events_{b}_before_{a}_after_snr_{s}'.format(b=before, a=after, s=suffix), tempcsat)
		Draw1DHistogram('strip_location_not_sat_events_{b}_before_{a}_after_snr_{s}'.format(b=before, a=after, s=suffix), tempcnosat)

	def DoSaturationAnalysis(self, cells='all', before=0, after=1):
		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz
		self.GetCutsFromCutManager(cells)
		self.GetVarzFromTranspGrid()
		self.DoProfileMaps(before=0, after=1)
		self.DoCellMaps(cells, before=0, after=1)
		self.PlotStripHistograms(cells, before=0, after=1)
		self.Do1DPHHistos(cells, before=0, after=1)
		if before != 0 or after != 1:
			self.PlotStripHistograms(cells, before=before, after=after)
			self.Do1DPHHistos(cells, before=before, after=after)
		# self.Do1DChSignalHistos(cells, False)
		# self.DoProfileMaps()
		# self.DoCellMaps(cells)
		# self.PlotStripHistogram(cells)

	def DefineRegion(self, before=0, after=1):
		if self.trans_grid.trans_tree.GetFriend('satRegions'):
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('satRegions'))
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=after, skipBefore=before)

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


		self.sat_adc_ch_cut = self.trans_grid.cuts_man.sat_adc_ch
		self.sat_adc_h_cut = self.trans_grid.cuts_man.sat_adc_h
		self.not_sat_adc_ch_cut = self.trans_grid.cuts_man.not_sat_adc_ch
		self.not_sat_adc_h_cut = self.trans_grid.cuts_man.not_sat_adc_h
		self.sat_adc_N_ch_cut = self.trans_grid.cuts_man.sat_adc_N_ch
		self.sat_adc_N_h_cut = self.trans_grid.cuts_man.sat_adc_N_h
		self.not_sat_adc_N_ch_cut = self.trans_grid.cuts_man.not_sat_adc_N_ch
		self.not_sat_adc_N_h_cut = self.trans_grid.cuts_man.not_sat_adc_N_h

		self.sat_evts_region = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region = self.trans_grid.cuts_man.not_sat_evts_region

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

if __name__ == '__main__':
	c = SaturationAnalysis(None, 0, 0)
