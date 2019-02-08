#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

tests = range(0, 10) + [100]
color_index = 10000

class TestAreas:
	def __init__(self, num=0, clust_size=1, dir='.', run=0, cellsize=50, do_fit=True):
		self.do_fit = do_fit
		self.num = num
		self.clust_size = clust_size
		self.dir = dir
		self.run = run
		self.cellsize = cellsize
		self.trans_grid = TransparentGrid(dir, run, cellsize)
		self.trans_grid.pkl_sbdir = 'test' + str(num)
		self.threshold = 800 if cellsize == 50 else 200
		self.window_shift = 10
		self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.5
		self.min_snr, self.max_snr = -64.25, 64.25
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.num_strips = self.trans_grid.num_strips if self.trans_grid.num_strips != 0 else 3
		self.cluster_size = self.trans_grid.num_strips if self.trans_grid.num_strips != 0 else 3

	def SetTest(self):
		self.trans_grid.cell_resolution = 50.0 / 13.0 if self.trans_grid.col_pitch == 50 else 100.0 / 51
		if self.num == 0:
			self.trans_grid.cell_resolution = 50.0 / 25.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(2 + 1 - self.clust_size, 17, 20)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 7, 9)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 17, 19)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 16, 20)
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 7, 9)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 17)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 8, 8)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 12, 12)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 21, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 24, 24)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 12, 12)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 17)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 19, 25)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 12, 17)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 19, 19)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 13, 17)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 19, 20)
			self.trans_grid.RemoveFromGoodArea(3 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 7)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 19)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 1:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			# self.trans_grid.cell_resolution = 50.0 / 17.0
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreas(2 + 1 - self.clust_size, 8)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 7, 9)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 16)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 17)
			# self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 17)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 17)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 12, 17)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 13, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 19)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 17)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 2:
			self.trans_grid.cell_resolution = 50.0 / 15.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 17)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 15, 19)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 3:
			self.trans_grid.cell_resolution = 50.0 / 15.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 4:
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 16)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 5:
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 19)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 18)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 6:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 14, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 14, 18)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 18)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 7:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 15, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 15, 17)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 8:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasRow(15, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddGoodAreasRow(16, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddGoodAreasRow(17, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 9:
			self.trans_grid.cell_resolution = 2
			self.trans_grid.ResetAreas()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterCharge{n}'.format(n=self.clust_size))
		elif self.num == 100:
			self.trans_grid.cell_resolution = 2
			self.trans_grid.ResetAreas()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterCharge{n}'.format(n=self.clust_size))

	def PlotTestClusterStudies(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_clusters = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]
		if cells == 'good':
			list_cuts_noise.append(self.trans_grid.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts_noise.append(self.trans_grid.gridAreas.badAreasCutNames_diamond)
		temp_cut_clusters = {i: '&&'.join(list_cut) for i, list_cut in list_cuts_clusters.iteritems()}
		temp_cut_noise = '&&'.join(list_cuts_noise)
		lastbin = int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
		tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
		suffix = 'not_selection' if cells == 'bad' else 'selection' if cells == 'good' else 'all'
		if cells == 'all':
			temp = [self.trans_grid.DrawPH('ph_c{n}_{c}'.format(n=i, c=suffix), self.min_snr, self.max_snr, self.delta_snr, var='diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='PH closestStrip{n} (SNR)'.format(n=i), cuts=temp_cut_clusters[i]) for i in xrange(self.cluster_size)]
		elif cells == 'good':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr, self.max_snr, int(np.floor((self.max_snr - self.min_snr) / float(self.delta_snr) + 0.5))
			temp = [self.trans_grid.DrawPHGoodAreas('ph_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		elif cells == 'bad':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr, self.max_snr, int(np.floor((self.max_snr - self.min_snr) / float(self.delta_snr) + 0.5))
			temp = [self.trans_grid.DrawPHBadAreas('ph_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins
		temp = [self.trans_grid.DrawPH('signal_noise_c{n}_{c}'.format(n=i, c=suffix), self.min_snr, self.max_snr, self.delta_snr, 'diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='goff') for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetLineColor(ro.kGray + 1) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].Scale(self.trans_grid.histo['ph_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)/self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetStats(0) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetTitle('signal_not_in_cluster(scaled)') for i in xrange(self.cluster_size)]
		for c in xrange(self.cluster_size):
			ro.gPad.Update()
			self.trans_grid.canvas['ph_c{n}_{c}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.canvas['ph_c{n}_{c}'.format(n=c, c=suffix)].cd()
			self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=c, c=suffix)].Draw('same')
			ro.gPad.Update()
			if self.trans_grid.histo['ph_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats'):
				self.trans_grid.histo['ph_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX1NDC(0.15)
				self.trans_grid.histo['ph_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX2NDC(0.45)
				ro.gPad.Update()
			legend = self.trans_grid.canvas['ph_c{n}_{c}'.format(n=c, c=suffix)].BuildLegend()
			ro.gPad.Update()
			legend.SetX1NDC(0.15)
			legend.SetX2NDC(0.45)
			legend.SetY1NDC(0.5)
			legend.SetY2NDC(0.6)
			ro.gPad.Update()
			self.w += self.window_shift

		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)] = ro.TCanvas('ph_ch{c}'.format(c=suffix), 'ph_ch{c}'.format(c=suffix), 1)
		for c in xrange(self.cluster_size):
			self.trans_grid.histo['ph_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetStats(0)
			self.trans_grid.histo['ph_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].Draw('same')
			r, g, b = ReturnRGB(c, 0, self.cluster_size)
			self.trash.append(ro.TColor(color_index + c, r, g, b))
			self.trans_grid.histo['ph_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetLineColor(color_index + c)
			self.trans_grid.histo['ph_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetMarkerColor(ro.kBlack)

		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].SetLogy()
		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].SetGridx()
		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].SetGridy()
		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].SetTicky()

		legend = self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].BuildLegend()
		ro.gPad.Update()
		legend.SetX1NDC(0.15)
		legend.SetX2NDC(0.45)
		legend.SetY1NDC(0.7)
		legend.SetY2NDC(0.9)
		ro.gPad.Update()
		self.trans_grid.canvas['ph_ch{c}'.format(c=suffix)].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

		for clch in ['1', '2', 'N']:
			for c in xrange(self.cluster_size):
				cuts = '({c})'.format(c=self.trans_grid.gridAreas.goodAreasCutNames_diamond) if cells == 'good' else '({c})'.format(c=self.trans_grid.gridAreas.badAreasCutNames_diamond) if cells == 'good' else ''
				self.trans_grid.DrawHisto2D('ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix), -500, 2500, 50, 'ph_ch{c}[ADC]'.format(c=c), 0, 4200, 40, 'ph{clch}[ADC]'.format(clch=clch), 'diaChSignal[clusterChannel{c}]'.format(c=c), 'clusterCharge{clch}'.format(clch=clch), cuts)
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetGridx()
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetGridy()
				ro.gPad.Update()
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift

	def PlotSaturation(self):
		self.trans_grid.DrawPHGoodAreas('PH_Saturated_Events', 'clusterCharge2', '((transparentEvent)&&((diaChADC[clusterChannel0]==4095)||(diaChADC[clusterChannel1]==4095)||(diaChADC[clusterChannel2]==4095)))')
		self.trans_grid.canvas['PH_Saturated_Events'].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

	def PlotTestForNegative(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_clusters = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr_neg, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr_neg)]
		if cells == 'good':
			list_cuts_noise.append(self.trans_grid.gridAreas.goodAreasCutNames_diamond)
		elif cells == 'bad':
			list_cuts_noise.append(self.trans_grid.gridAreas.badAreasCutNames_diamond)
		temp_cut_clusters = {i: '&&'.join(list_cut) for i, list_cut in list_cuts_clusters.iteritems()}
		temp_cut_noise = '&&'.join(list_cuts_noise)
		lastbin = int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
		tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
		suffix = 'not_selection' if cells == 'bad' else 'selection' if cells == 'good' else 'all'
		if cells == 'all':
			temp = [self.trans_grid.DrawPH('ph_neg_c{n}_{c}'.format(n=i, c=suffix), self.min_snr_neg, self.max_snr_neg, self.delta_snr, var='diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='PH closestStrip{n} (SNR)'.format(n=i), cuts=temp_cut_clusters[i]) for i in xrange(self.cluster_size)]
		elif cells == 'good':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr_neg, self.max_snr_neg, int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
			temp = [self.trans_grid.DrawPHGoodAreas('ph_neg_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		elif cells == 'bad':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr_neg, self.max_snr_neg, int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
			temp = [self.trans_grid.DrawPHBadAreas('ph_neg_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins
		temp = [self.trans_grid.DrawPH('signal_noise_c{n}_{c}'.format(n=i, c=suffix), self.min_snr_neg, self.max_snr_neg, self.delta_snr, 'diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='goff') for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetLineColor(ro.kGray + 1) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].Scale(self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)/self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetStats(0) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetTitle('signal_not_in_cluster(scaled)') for i in xrange(self.cluster_size)]
		for c in xrange(self.cluster_size):
			ro.gPad.Update()
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].cd()
			self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=c, c=suffix)].Draw('same')
			ro.gPad.Update()
			if self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats'):
				self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX1NDC(0.15)
				self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX2NDC(0.45)
				ro.gPad.Update()
			hbla = ro.TH1F('bla{n}_{c}'.format(n=c, c=suffix), 'Negative cut', 2, self.max_snr_neg + 1, self.max_snr_neg + 2)
			hbla.SetLineStyle(7)
			hbla.SetLineColor(ro.kRed)
			hbla.SetLineWidth(2)
			hbla.Draw('same')
			self.trash.append(hbla)
			ro.gPad.Update()
			legend = self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].BuildLegend()
			ro.gPad.Update()
			legend.SetX1NDC(0.15)
			legend.SetX2NDC(0.45)
			legend.SetY1NDC(0.5)
			legend.SetY2NDC(0.6)
			ro.gPad.Update()
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)] = ro.TLine(-self.trans_grid.neg_cut, 0, -self.trans_grid.neg_cut, self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].GetMaximum())
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineColor(ro.kRed)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineStyle(7)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineWidth(2)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].Draw('same')
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLogy()
			ro.gPad.Update()
			self.w += self.window_shift

		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)] = ro.TCanvas('ph_neg_{c}'.format(c=suffix), 'ph_neg_{c}'.format(c=suffix), 1)
		for c in xrange(self.cluster_size):
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetStats(0)
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].Draw('same')
			r, g, b = ReturnRGB(c, 0, self.cluster_size)
			self.trash.append(ro.TColor(color_index + c, r, g, b))
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetLineColor(color_index + c)
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetMarkerColor(ro.kBlack)

		hbla = ro.TH1F('bla000', 'Negative cut', 2, self.max_snr_neg + 1, self.max_snr_neg + 2)
		hbla.SetLineStyle(7)
		hbla.SetLineColor(ro.kRed)
		hbla.SetLineWidth(2)
		hbla.Draw('same')
		self.trash.append(hbla)
		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].SetLogy()
		legend = self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].BuildLegend()
		ro.gPad.Update()
		legend.SetX1NDC(0.15)
		legend.SetX2NDC(0.45)
		legend.SetY1NDC(0.7)
		legend.SetY2NDC(0.9)
		ro.gPad.Update()
		self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=self.cluster_size - 1, c=suffix)].Draw('same')
		ro.gPad.Update()
		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

		for c in xrange(self.cluster_size):
			self.trans_grid.DrawProfile2DDiamond('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diaChSignal', '(diaChannels==clusterChannel{i})'.format(i=c))
			self.trans_grid.profile['ph_neg_map_c{n}'.format(n=c, c=suffix)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph_neg_map_c{n}'.format(n=c, c=suffix)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph_neg_map_c{n}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawTCutGs('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters('ph_neg_map_c{n}'.format(n=c, c=suffix))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DNoTopBottomBorders('signal_neg_map_c{n}'.format(n=c, c=suffix), -0.5, 127.5, 1, 'dia X ch', self.trans_grid.row_info_diamond['0']-np.floor(self.trans_grid.row_info_diamond['0'] / self.trans_grid.row_info_diamond['pitch'] + 0.5) * self.trans_grid.row_info_diamond['pitch'],
			                                                self.trans_grid.row_info_diamond['0'] + (256 - np.floor(self.trans_grid.row_info_diamond['0'] / self.trans_grid.row_info_diamond['pitch'] + 0.5)) * self.trans_grid.row_info_diamond['pitch'], float(self.trans_grid.row_info_diamond['pitch'])/self.trans_grid.bins_per_ch_y,
			                                                'sil pred Y [#mum]', 'diaChannels', 'diaChYPred', 'diaChSignal', 'Ch Signal [ADC]', '(diaChannels==clusterChannel{i})&&(diaChSignal/diaChPedSigmaCmc<=-{c})'.format(i=c, c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['signal_neg_map_c{n}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.profile['signal_neg_map_c{n}'.format(n=c, c=suffix)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['signal_neg_map_c{n}'.format(n=c, c=suffix)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs('signal_neg_map_c{n}'.format(n=c, c=suffix), 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters('signal_neg_map_c{n}'.format(n=c, c=suffix))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph_neg_cells_c{n}_{c}'.format(n=c, c=suffix), 'diaChSignal', cells, '(diaChannels==clusterChannel{n})&&(diaChSignal/diaChPedSigmaCmc<=-{c})'.format(n=c, c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph_neg_cells_c{n}_{c}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift


	def PlotTest(self):
		num = self.num
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		tempn = self.num_strips if self.num_strips == 1 else 2
		for ch in xrange(1, tempn + 1):
			#  plot map
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts='({l}<=diaChYPred)&&(diaChYPred<={h})'.format(l=y0, h=yup))
			# ro.gPad.Update()
			self.trans_grid.canvas['ph{c}_map_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.gridTextDiamond.Draw('same TEXT0')
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			#  draw selected areas
			self.trans_grid.DrawGoodAreasDiamond('ph{c}_map_test{n}'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}'.format(c=ch, n=num))
			#  plot cell overlay
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good')
			self.trans_grid.canvas['ph{c}_cells_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}'.format(c=ch, n=num), 'goff')
			#  show center region in cell
			#  Efficiency plots
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts='(clusterCharge{c}>=5*diaChPedSigmaCmc[clusterChannel0])'.format(c=ch), plot_option='prof colz goff')
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)] = self.trans_grid.histo['hit_map_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Clone('h_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetTitle('h_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Sumw2()
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Divide(self.trans_grid.histo['hit_map_ph{c}_cells_test{n}'.format(c=ch, n=num)])
			self.trans_grid.canvas['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)] = ro.TCanvas('c_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), 'c_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), 1)
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Draw('colz')
			SetDefault2DStats(self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)])
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].GetZaxis().SetTitle('Efficiency')
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].GetZaxis().SetRangeUser(0.9, 1)
			self.trans_grid.canvas['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawEfficiencyADCCut('Eff_ph{c}VsADC_test{n}'.format(c=ch, n=num), 'clusterCharge' + str(ch), cells='good', cut='transparentEvent', ymin_plot=0.95)
			self.trans_grid.canvas['Eff_ph{c}VsADC_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			#  draw ph of selected areas
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch))
			self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=ch, n=num), cells='good', var='clusterCharge' + str(ch))
			self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			#  fit distribution for central region
			if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c=ch, n=num), color=ro.kRed)
			#  get difference between cell and center
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=ch, n=num))
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Reset('ICES')
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c=ch, n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c=ch, n=num)], 1, -1)
			self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			if self.do_fit:
				self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlue)
				self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].cd()
				self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c=ch, n=num)].fit.Draw('same')
				self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c=ch, n=num)].fit.Draw('same')
				self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c=ch, n=num), 'ph{c}_test{n}_centers'.format(c=ch, n=num), 'ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlack)
			# ro.gPad.Update()
			#  position of negative clusters
			cut_no_neg = '(Sum$((diaChHits)&&(diaChSignal>-{c}*diaChPedSigmaCmc))=={n})'.format(c=self.trans_grid.neg_cut, n=self.cluster_size)
			cut_any_neg = '(Sum$((diaChHits)&&(diaChSignal<-{c}*diaChPedSigmaCmc))>0)'.format(c=self.trans_grid.neg_cut)
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts=cut_any_neg)
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts=cut_any_neg)
			self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cuts=cut_any_neg)
			self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=ch, n=num), color=ro.kBlue)
			self.w += self.window_shift

			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts=cut_no_neg)
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts=cut_no_neg)
			self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cuts=cut_no_neg)
			self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=ch, n=num), color=ro.kBlue)
			self.w += self.window_shift

		if self.clust_size >= 2:
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1', cells='good')
			self.trans_grid.canvas['ph2_minus_ph1_map_test{n}'.format(n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawTCutCentersInCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num))
			tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.trans_grid.phmin_neg, self.trans_grid.phmax_neg, self.trans_grid.phbins_neg
			self.trans_grid.DrawPHGoodAreas('ph2_minus_ph1_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1')
			self.trans_grid.canvas['ph2_minus_ph1_test{n}'.format(n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins

			if self.num_strips != 2:
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num), 'clusterChargeN-clusterCharge2', cells='good')
				self.trans_grid.canvas['ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawTCutCentersInCellOverlay('ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num))
				tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
				self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.trans_grid.phmin_neg, self.trans_grid.phmax_neg, self.trans_grid.phbins_neg
				self.trans_grid.DrawPHGoodAreas('ph{s}_minus_ph2_test{n}'.format(s=self.num_strips, n=num), 'clusterChargeN-clusterCharge2')
				self.trans_grid.canvas['ph{s}_minus_ph2_test{n}'.format(s=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins

				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good')
				self.trans_grid.canvas['ph{c}_cells_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				#  show center region in cell
				self.trans_grid.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=self.num_strips, n=num))
				#  draw ph of selected areas
				self.trans_grid.DrawEfficiencyADCCut('Eff_ph{c}VsADC_test{n}'.format(c=self.num_strips, n=num), 'clusterChargeN', cells='good', cut='transparentEvent', ymin_plot=0.95)
				self.trans_grid.canvas['Eff_ph{c}VsADC_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c=self.num_strips, n=num), var='clusterChargeN')
				self.trans_grid.canvas['ph{c}_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), cells='good', var='clusterChargeN')
				self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				#  fit distribution for central region
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), color=ro.kRed)
				#  get difference between cell and center
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num))
				self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].Reset('ICES')
				self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c=self.num_strips, n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)], 1, -1)
				self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				if self.do_fit:
					self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num), color=ro.kBlue)
					self.trans_grid.canvas['ph{c}_test{n}'.format(c=self.num_strips, n=num)].cd()
					self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)].fit.Draw('same')
					self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].fit.Draw('same')
					self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c=self.num_strips, n=num), 'ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), 'ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num), color=ro.kBlack)
				# ro.gPad.Update()
				#  position of negative clusters
				cut_no_neg = '(Sum$((diaChHits)&&(diaChSignal>-{c}*diaChPedSigmaCmc))=={n})'.format(c=self.trans_grid.neg_cut, n=self.cluster_size)
				cut_any_neg = '(Sum$((diaChHits)&&(diaChSignal<-{c}*diaChPedSigmaCmc))>0)'.format(c=self.trans_grid.neg_cut)
				self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num), varz='clusterChargeN', cuts=cut_any_neg)
				self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.w += self.window_shift
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good', cuts=cut_any_neg)
				self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cuts=cut_any_neg)
				self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=self.num_strips, n=num), color=ro.kBlue)
				self.w += self.window_shift

				self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num), varz='clusterChargeN', cuts=cut_no_neg)
				self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.w += self.window_shift
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good', cuts=cut_no_neg)
				self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cuts=cut_no_neg)
				self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num), color=ro.kBlue)
				self.w += self.window_shift

	def SaveCanvas(self):
		self.trans_grid.SaveCanvasInlist(self.trans_grid.canvas.keys())

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
	parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
	parser.add_option('-c', '--cellsize', dest='cellsize', type='int', default=50, help='cell size of the square 3D device')
	parser.add_option('-t', '--test', dest='testnumber', type='int', default=-1, help='Run a automatically one of the predefined tests')
	parser.add_option('-f', '--dofit', dest='dofit', default=False, action='store_true', help='Enables fitting')
	parser.add_option('-n', '--numstrips', dest='numstrips', type='int', default=2, help='Number of strips to use')
	parser.add_option('-a', '--auto', dest='auto', default=False, action='store_true', help='Sets up test, creates plots and saves them automatically if toggled')
	parser.add_option('-s', '--saturation', dest='saturation', default=False, action='store_true', help='Sets up saturation plots')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	testnum = int(options.testnumber)
	do_fit = bool(options.dofit)
	numstrips = int(options.numstrips)
	cellsize = int(options.cellsize) if testnum < 100 else 100
	do_sat = bool(options.saturation)
	autom = bool(options.auto)

	t = TestAreas(testnum, numstrips, dir, run, cellsize, do_fit)
	t.trans_grid.SetLines()
	t.trans_grid.CreateTCutGs()
	if testnum in tests:
		t.SetTest()
		if t.trans_grid.loaded_pickle:
			t.trans_grid.LoadPickle()
		else:
			t.trans_grid.FindPickleValues()
			ExitMessage('Run it again to load the generated pickle :)', os.EX_OK)
		if t.trans_grid.loaded_default_pickle:
			t.trans_grid.FindBinningAndResolution()
			t.trans_grid.SavePickle()
		if autom:
			t.PlotTestClusterStudies('good')
			t.PlotTestForNegative('good')
			t.PlotTest()
			if do_sat:
				t.PlotSaturation()
			t.SaveCanvas()

