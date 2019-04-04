#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class NegativeChargesAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize):
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

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

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

	def OverlayNoiseDistribution(self, histo, cells='all'):
		suffix = self.suffix[cells]
		hname = histo.GetName().split('h_')[1]
		typ = 'adc' if 'adc' in hname.lower() else 'snr'
		noise_name0 = 'signal_noise_{s}_{t}'.format(s=suffix, t=typ)
		if not noise_name0 in self.trans_grid.histo.keys():
			self.PlotNoiseNotInCluster(cells)
		elif not self.trans_grid.histo[noise_name0]:
			del self.trans_grid.histo[noise_name0]
			self.PlotNoiseNotInCluster(cells)

		noise_name_new = noise_name0 + '_' + hname
		nbins = histo.GetNbinsX()
		xbins = np.zeros(nbins, 'float64')
		histo.GetXaxis().GetLowEdge(xbins)
		xbins = np.append(xbins, 2 * xbins[-1] - xbins[-2])
		self.trans_grid.histo[noise_name_new] = self.trans_grid.histo[noise_name0].Rebin(nbins, 'h_' + noise_name_new, xbins)
		if self.trans_grid.histo[noise_name_new]:
			self.trans_grid.histo[noise_name_new].SetTitle(noise_name0 + '(scaled)')
			scale = histo.GetMaximum() / self.trans_grid.histo[noise_name_new].GetMaximum()
			self.trans_grid.histo[noise_name_new].Scale(scale)
			self.trans_grid.histo[noise_name_new].SetLineColor(ro.kGray + 1)
			self.trans_grid.histo[noise_name_new].SetStats(0)
			self.trans_grid.canvas[hname].cd()
			if self.trans_grid.histo[noise_name_new].GetFunction('f_gaus_' + noise_name0):
				self.trans_grid.histo[noise_name_new].GetFunction('f_gaus_' + noise_name0).SetBit(ro.TF1.kNotDraw)
			self.trans_grid.histo[noise_name_new].Draw('same')
			ro.gPad.Update()

	def Do1DChSignalHistos(self, cells='all', doLog=False):
		def DrawHisto(name, histo_limits, plot_lims, deltax, varz, varname, cuts):
			self.trans_grid.DrawHisto1D(name, histo_limits['min'], histo_limits['max'], deltax, varz, varname, cuts)
			self.trans_grid.histo[name].GetXaxis().SetRangeUser(plot_lims['min'], plot_lims['max'])
			SetX1X2NDC(self.trans_grid.histo[name], 0.15, 0.45, 'stats')
			self.OverlayNoiseDistribution(self.trans_grid.histo[name], cells)
			if doLog: self.trans_grid.canvas[name].SetLogy()
			legend = self.trans_grid.canvas[name].BuildLegend()
			ro.gPad.Update()
			SetLegendX1X2Y1Y2(legend, 0.15, 0.45, 0.5, 0.6)
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		suffix = suffix + '_logScale' if doLog else suffix

		tempcutsadc = self.neg_adc_phN_h['PH{c}_H'.format(c=self.cluster_size)]
		tempcutssnr = self.neg_snr_phN_h['PH{c}_H'.format(c=self.cluster_size)]
		for ch in xrange(self.cluster_size):
			if 'PH_Ch' + str(ch) in self.ph_snr_ch_varz.keys():
				tempcuts = tempcutssnr
				minz, maxz = self.min_snr_neg, self.max_snr_neg
				hist_limits = Get1DLimits(minz, maxz, self.delta_snr, 2)
				plot_limits = Get1DLimits(minz, maxz, self.delta_snr, 1)
				DrawHisto('PH_Ch{c}_snr_neg_evts_{s}'.format(c=ch, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_ch_varz['PH_Ch{c}'.format(c=ch)], 'PH cluster ch{c} neg events [SNR]'.format(c=ch), tempcuts)

			if 'PH_Ch' + str(ch) in self.ph_adc_ch_varz.keys():
				tempcuts = tempcutsadc
				minz, maxz = self.min_adc_neg, self.max_adc_neg
				hist_limits = Get1DLimits(minz, maxz, self.delta_adc, 2)
				plot_limits = Get1DLimits(minz, maxz, self.delta_adc, 1)
				DrawHisto('PH_Ch{c}_adc_neg_evts_{s}'.format(c=ch, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_ch_varz['PH_Ch{c}'.format(c=ch)], 'PH cluster ch{c} neg events [ADC]'.format(c=ch), tempcuts)

			if 'PH_H' + str(ch+1) in self.ph_snr_h_varz.keys():
				tempcuts = tempcutssnr
				minz, maxz = self.min_snr_neg, self.max_snr_neg
				hist_limits = Get1DLimits(minz, maxz, self.delta_snr, 2)
				plot_limits = Get1DLimits(minz, maxz, self.delta_snr, 1)
				DrawHisto('PH_H{c}_snr_neg_evts_{s}'.format(c=ch+1, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_h_varz['PH_H{c}'.format(c=ch+1)], 'PH highest {c}{sf} ch neg events [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

			if 'PH_H' + str(ch+1) in self.ph_adc_h_varz.keys():
				tempcuts = tempcutsadc
				minz, maxz = self.min_adc_neg, self.max_adc_neg
				hist_limits = Get1DLimits(minz, maxz, self.delta_adc, 2)
				plot_limits = Get1DLimits(minz, maxz, self.delta_adc, 1)
				DrawHisto('PH_H{c}_adc_neg_evts_{s}'.format(c=ch+1, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_h_varz['PH_H{c}'.format(c=ch+1)], 'PH highest {c}{sf} ch neg events [ADC]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

	def DoProfileMaps(self):
		xmin, xmax, deltax, xname = self.trans_grid.ch_ini - 1.5, self.trans_grid.ch_end + 1.5, 1.0/self.trans_grid.bins_per_ch_x, 'pred dia hit ch',
		ymin, ymax, deltay, yname = self.trans_grid.row_info_diamond['0'] - RoundInt(float(self.trans_grid.row_info_diamond['0']) / self.trans_grid.row_info_diamond['pitch'], 'f8') * self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['0'] + (256 - RoundInt(float(self.trans_grid.row_info_diamond['0']) / self.trans_grid.row_info_diamond['pitch'], 'f8')) * self.trans_grid.row_info_diamond['pitch'], float(self.trans_grid.row_info_diamond['pitch'])/self.trans_grid.bins_per_ch_y, 'sil pred dia hit in Y [#mum]'

		def DrawProfile2D(name, varz, varzname, cut, xdelt=deltax, namex=xname, varx='diaChXPred'):
			self.trans_grid.DrawProfile2D(name, xmin, xmax, xdelt, namex, ymin, ymax, deltay, yname, varx, 'diaChYPred', varz, varzname, cut)
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)

		for ch in xrange(self.cluster_size):
			tempc = self.neg_snr_ph_ch['PH_Ch' + str(ch)]
			DrawProfile2D('PH_Ch{c}_neg_map_pred_hit_snr'.format(c=ch), self.ph_snr_ch_varz['PH_Ch' + str(ch)], 'PH cluster ch{c} [SNR]'.format(c=ch), tempc)
			tempc = self.neg_adc_ph_ch['PH_Ch' + str(ch)]
			DrawProfile2D('PH_Ch{c}_neg_map_pred_hit_adc'.format(c=ch), self.ph_adc_ch_varz['PH_Ch' + str(ch)], 'PH cluster ch{c} [ADC]'.format(c=ch), tempc)

			tempc = self.neg_snr_ph_h['PH_H' + str(ch+1)]
			DrawProfile2D('PH_H{c}_neg_map_pred_hit_snr'.format(c=ch+1), self.ph_snr_h_varz['PH_H' + str(ch+1)], 'PH highest {c}{sf} ch [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc)
			tempc = self.neg_adc_ph_h['PH_H' + str(ch+1)]
			DrawProfile2D('PH_H{c}_neg_map_pred_hit_adc'.format(c=ch+1), self.ph_adc_h_varz['PH_H' + str(ch+1)], 'PH highest {c}{sf} ch [ADC]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc)

			tempc = self.neg_snr_ph_ch['PH_Ch' + str(ch)]
			DrawProfile2D('PH_Ch{c}_neg_map_pred_ch_snr'.format(c=ch), self.ph_snr_ch_varz['PH_Ch' + str(ch)], 'PH cluster ch{c} [SNR]'.format(c=ch), tempc, 1, 'cluster ch{c}'.format(c=ch), 'clusterChannel{c}'.format(c=ch))
			tempc = self.neg_adc_ph_ch['PH_Ch' + str(ch)]
			DrawProfile2D('PH_Ch{c}_neg_map_pred_ch_adc'.format(c=ch), self.ph_adc_ch_varz['PH_Ch' + str(ch)], 'PH cluster ch{c} [ADC]'.format(c=ch), tempc, 1, 'cluster ch{c}'.format(c=ch), 'clusterChannel{c}'.format(c=ch))

			tempc = self.neg_snr_ph_h['PH_H' + str(ch+1)]
			DrawProfile2D('PH_H{c}_neg_map_pred_ch_snr'.format(c=ch+1), self.ph_snr_h_varz['PH_H' + str(ch+1)], 'PH highest {c}{sf} ch [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc, 1, 'highest {c}{sf} ch'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 'clusterChannelHighest{c}'.format(c=ch+1))
			tempc = self.neg_adc_ph_h['PH_H' + str(ch+1)]
			DrawProfile2D('PH_H{c}_neg_map_pred_ch_adc'.format(c=ch+1), self.ph_adc_h_varz['PH_H' + str(ch+1)], 'PH highest {c}{sf} ch [ADC]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc, 1, 'highest {c}{sf} ch'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 'clusterChannelHighest{c}'.format(c=ch+1))

	def Do1DPHHistos(self, cells='all'):
		suffix = self.suffix[cells]
		noise_name0 = 'signal_noise_{s}_{t}'.format(s=suffix, t='adc')
		sigm = self.trans_grid.histo[noise_name0].GetRMS()
		def DrawPHHisto(name, varz, varzname, cuts, xmin=10000, xmax=-10000, deltax=-1):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varzname, xmin=xmin, xmax=xmax, deltax=deltax)
			self.PosCanvas(name)

		tempcsnr = self.neg_adc_phN_h['PH{c}_H'.format(c=self.cluster_size)]
		tempcadc = self.neg_adc_phN_h['PH{c}_H'.format(c=self.cluster_size)]
		minsnr, maxsnr = int(RoundInt(self.trans_grid.phmin / sigm)), int(RoundInt(self.trans_grid.phmax / sigm))
		deltsnr = float(maxsnr - minsnr) / float(self.trans_grid.phbins)
		for ch in xrange(1, self.cluster_size + 1):
			tempc = tempcsnr
			DrawPHHisto('PH{c}_Ch_neg_events_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} cluster chs [SNR]', tempc, minsnr, maxsnr, deltsnr)
			tempc = tempcadc
			DrawPHHisto('PH{c}_Ch_neg_events_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} cluster chs [ADC]', tempc)

			tempc = tempcsnr
			DrawPHHisto('PH{c}_H_neg_events_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [SNR]', tempc, minsnr, maxsnr, deltsnr)
			tempc = tempcadc
			DrawPHHisto('PH{c}_H_neg_events_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [ADC]', tempc)

	def DoNegativeAnalysis(self, cells='all'):
		self.GetCutsFromCutManager(cells)
		self.GetVarzFromTranspGrid()
		self.Do1DChSignalHistos(cells, False)
		self.DoProfileMaps()
		self.Do1DPHHistos(cells)
		# self.DoProfileMaps()
		# self.DoPedestalEventHistograms(False)
		# self.DoStrips2DHistograms()
		# self.PlotNoiseNotInCluster(cells)
		# self.PlotNoiseNCChannels('all')

	def GetCutsFromCutManager(self, cells):
		self.noise_cuts[cells] = self.trans_grid.cuts_man.noise_cuts[cells]
		self.noise_friend_cuts[cells] = self.trans_grid.cuts_man.noise_friend_cuts[cells]
		self.in_transp_cluster = self.trans_grid.cuts_man.ConcatenateCuts(cut2=self.trans_grid.cuts_man.in_transp_cluster, cut1=self.trans_grid.cuts_man.transp_ev)
		self.noise_nc_cuts[cells] = self.trans_grid.cuts_man.noise_nc_cuts[cells]
		self.noise_nc_friend_cuts[cells] = self.trans_grid.cuts_man.noise_nc_friend_cuts[cells]

		self.neg_snr_ph_ch = self.trans_grid.cuts_man.neg_snr_ph_ch
		self.neg_snr_ph_h = self.trans_grid.cuts_man.neg_snr_ph_h
		self.not_neg_snr_ph_ch = self.trans_grid.cuts_man.not_neg_snr_ph_ch
		self.not_neg_snr_ph_h = self.trans_grid.cuts_man.not_neg_snr_ph_h

		self.neg_adc_ph_ch = self.trans_grid.cuts_man.neg_adc_ph_ch
		self.neg_adc_ph_h = self.trans_grid.cuts_man.neg_adc_ph_h
		self.not_neg_adc_ph_ch = self.trans_grid.cuts_man.not_neg_adc_ph_ch
		self.not_neg_adc_ph_h = self.trans_grid.cuts_man.not_neg_adc_ph_h

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

if __name__ == '__main__':
	c = NegativeChargesAnalysis(None, 0, 0)
