#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class TestAreas:
	def __init__(self, configfile='', run=0):
		self.do_fit = False
		self.do_saturation = True
		self.num = 0
		self.clust_size = 2
		self.dir = '.'
		self.run = run
		self.cellsize = 50
		self.threshold = 800
		self.do_threshold = False
		self.window_shift = 5
		self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.125
		self.min_snr, self.max_snr = -64.25, 64.25
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.num_rows = 0
		self.rows = []
		self.cols = []
		self.cells = []
		self.rcells = []
		self.config_file = configfile
		if self.config_file != '':
			self.config_file = Correct_Path(self.config_file)
			self.ReadConfigFile()
		self.trans_grid = TransparentGrid(self.dir, self.run, self.cellsize)
		self.trans_grid.pkl_sbdir = 'test' + str(self.num)
		self.bias = self.trans_grid.bias
		self.saturated_ADC = self.trans_grid.saturated_ADC
		self.num_strips = self.trans_grid.num_strips if self.trans_grid.num_strips != 0 else 3
		self.cluster_size = self.trans_grid.num_strips if self.trans_grid.num_strips != 0 else 3
		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}
		if self.num_rows != 0:
			self.trans_grid.row_info_diamond['num'] = self.num_rows

	def ReadConfigFile(self):
		def unpack_row_col(string):
			elements = string.replace('{', '').replace('}', '')
			elements = elements.split(';')
			elements = [elem.split(',') for elem in elements]
			elements = [[int(elemij) for elemij in elemi] for elemi in elements]
			return elements
		if os.path.isfile(self.config_file):
			pars = ConfigParser()
			pars.read(self.config_file)
			print 'Reading config file for test area...', ; sys.stdout.flush()
			if pars.has_section('SETTINGS'):
				if pars.has_option('SETTINGS', 'run'):
					self.run = pars.getint('SETTINGS', 'run')
				if pars.has_option('SETTINGS', 'dir'):
					self.dir = Correct_Path(pars.get('SETTINGS', 'dir'))
				if pars.has_option('SETTINGS', 'cluster_size'):
					self.cluster_size = pars.getint('SETTINGS', 'cluster_size')
				if pars.has_option('SETTINGS', 'cell_size'):
					self.cellsize = pars.getint('SETTINGS', 'cell_size')
				if pars.has_option('SETTINGS', 'do_fit'):
					self.do_fit = pars.getboolean('SETTINGS', 'do_fit')
				if pars.has_option('SETTINGS', 'do_saturation'):
					self.do_saturation = pars.getboolean('SETTINGS', 'do_saturation')
				if pars.has_option('SETTINGS', 'test_number'):
					self.num = pars.getint('SETTINGS', 'test_number')
				if pars.has_option('SETTINGS', 'threshold'):
					self.threshold = pars.getfloat('SETTINGS', 'threshold')
				if pars.has_option('SETTINGS', 'do_threshold'):
					self.do_threshold = pars.getboolean('SETTINGS', 'do_threshold')
			if pars.has_section('ROWS'):
				if pars.has_option('ROWS', 'rows'):
					rows = pars.get('ROWS', 'rows')
					self.rows = unpack_row_col(rows)
				if pars.has_option('ROWS', 'num'):
					self.num_rows = pars.getint('ROWS', 'num')
			if pars.has_section('COLUMNS'):
				if pars.has_option('COLUMNS', 'cols'):
					cols = pars.get('COLUMNS', 'cols')
					self.cols = unpack_row_col(cols)
			if pars.has_section('CELLS'):
				if pars.has_option('CELLS', 'cells'):
					cells = pars.get('CELLS', 'cells')
					self.cells = unpack_row_col(cells)
			if pars.has_section('REMOVECELLS'):
				if pars.has_option('REMOVECELLS', 'rcells'):
					rcells = pars.get('REMOVECELLS', 'rcells')
					self.rcells = unpack_row_col(rcells)
			print 'Done'

	def SetTest(self):
		if self.do_threshold:
			self.trans_grid.ResetAreas()
			print 'Selecting areas with a ph2 greater or equal than', self.threshold, '...', ; sys.stdout.flush()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterCharge2')
			print 'Done'
			self.trans_grid.AddRemainingToBadAreas()
			print 'Marked the remaining cells as bad'
			self.trans_grid.gridAreas.SimplifyGoodAndBadAreas()
		elif len(self.rows) + len(self.cols) + len(self.cells) > 0:
			self.trans_grid.ResetAreas()
			for row in self.rows:
				self.trans_grid.AddGoodAreasRow(row[0], row[1], row[2])
				print 'Added row {r} from {coli} to {colj} to selection'.format(r=row[0], coli=row[1], colj= row[2])
			for col in self.cols:
				self.trans_grid.AddGoodAreasCol(col[0], col[1], col[2])
				print 'Added column {r} from {coli} to {colj} to selection'.format(r=col[0], coli=col[1], colj=col[2])
			for cell in self.cells:
				self.trans_grid.AddGoodAreas(cell[0], cell[1])
				print 'Added cell with column {c} and row {r} to selection'.format(c=cell[0], r=cell[1])
			for rcell in self.rcells:
				self.trans_grid.RemoveFromGoodArea(rcell[0], rcell[1])
				print 'Removed cell with column {c} and row {r} from selection'.format(c=rcell[0], r=rcell[1])
			self.trans_grid.AddRemainingToBadAreas()
			self.trans_grid.gridAreas.SimplifyGoodAndBadAreas()
			print 'Marked the remaining cells as bad'
		else:
			print 'Enter a correct settings file for the test area in variable config_file and re run ReadConfigFile before setting the test...'

	def PlotNoiseNotInCluster(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_noise_snr = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(transparentEvent)&&(diaChHits==0)&&(diaChSeed==0)&&(diaChsScreened==0)&&(diaChsNoisy==0)&&(diaChsNC==0)&&(TMath::Abs(diaChSignal)/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]
		list_cuts_noise_adc = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(transparentEvent)&&(diaChHits==0)&&(diaChSeed==0)&&(diaChsScreened==0)&&(diaChsNoisy==0)&&(diaChsNC==0)&&(TMath::Abs(diaChSignal)<{m}*diaChPedSigmaCmc)&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]
		if cells == 'good':
			list_cuts_noise_snr.append(self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond)
		elif cells == 'bad':
			list_cuts_noise_snr.append(self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
		temp_cut_noise_snr = '&&'.join(list_cuts_noise_snr)
		temp_cut_noise_adc = '&&'.join(list_cuts_noise_adc)
		lastbin = int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
		tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
		temph = ro.TH1F('temph0', 'temph0', 100, -32, 32)
		self.trans_grid.trans_tree.Draw('diaChSignal>>temph0', temp_cut_noise_adc, 'goff')
		mean, sigma = temph.GetMean(), temph.GetRMS()
		temph.Delete()
		suffix = self.suffix[cells]
		self.trans_grid.DrawPH('signal_noise_{c}_snr'.format(c=suffix), self.min_snr, self.max_snr, self.delta_snr, 'diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='Signal not in cluster (SNR)', cuts=temp_cut_noise_snr, option='goff')
		self.trans_grid.DrawPH('signal_noise_{c}_adc'.format(c=suffix), self.min_snr * sigma, self.max_snr * sigma, self.delta_snr * sigma, 'diaChSignal', varname='Signal not in cluster (ADC)', cuts=temp_cut_noise_snr, option='goff')

	def PlotTestClusterStudies(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_clusters_snr_ci = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise_snr_ci = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]
		list_cuts_clusters_adc_ci = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal<{m}*diaChPedSigmaCmc)&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise_adc_ci = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal<{m}*diaChPedSigmaCmc)&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]
		list_cuts_clusters_adc_phjk = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal<{m}*diaChPedSigmaCmc)&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise_adc_phjk = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal<{m}*diaChPedSigmaCmc)&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr)]

		if cells == 'good':
			list_cuts_noise_snr_ci.append(self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond)
		elif cells == 'bad':
			list_cuts_noise_snr_ci.append(self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
		temp_cut_clusters = {i: '&&'.join(list_cut) for i, list_cut in list_cuts_clusters_snr_ci.iteritems()}
		temp_cut_noise = '&&'.join(list_cuts_noise_snr_ci)
		lastbin = int(np.floor((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr) + 0.5))
		tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
		suffix = self.suffix[cells]
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
				cuts = '({c})'.format(c=self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond) if cells == 'good' else '({c})'.format(c=self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond) if cells == 'good' else ''
				self.trans_grid.DrawHisto2D('ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix), -500, 2500, 50, 'ph_ch{c}[ADC]'.format(c=c), 0, 4200, 40, 'ph{clch}[ADC]'.format(clch=clch), 'diaChSignal[clusterChannel{c}]'.format(c=c), 'clusterCharge{clch}'.format(clch=clch), cuts)
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetGridx()
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetGridy()
				ro.gPad.Update()
				self.trans_grid.canvas['ph_ch{ch}_vs_ph{clch}_{c}'.format(ch=c, clch=clch, c=suffix)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift

		for c in xrange(self.cluster_size):
			self.trans_grid.DrawProfile2DDiamond('ph_map_c{n}'.format(n=c, c=suffix), 'diaChSignal[clusterChannel{i}]'.format(i=c))
			self.trans_grid.profile['ph_map_c{n}'.format(n=c, c=suffix)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph_map_c{n}'.format(n=c, c=suffix)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph_map_c{n}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawTCutGs('ph_map_c{n}'.format(n=c, c=suffix), 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters('ph_map_c{n}'.format(n=c, c=suffix))
			self.w += self.window_shift

	def PlotSaturation(self):
		for c in xrange(1, self.cluster_size):
			self.trans_grid.DrawPHGoodAreas('ph{c}_saturated'.format(c=c), 'clusterCharge{c}'.format(c=c), '((transparentEvent)&&((diaChADC[clusterChannel0]=={s})||(diaChADC[clusterChannel1]=={s})||(diaChADC[clusterChannel2]=={s})))'.format(s=self.saturated_ADC))
			self.trans_grid.canvas['ph{c}_saturated'.format(c=c)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_not_saturated'.format(c=c), 'clusterCharge{c}'.format(c=c), '((transparentEvent)&&((diaChADC[clusterChannel0]<{s})&&(diaChADC[clusterChannel1]<{s})&&(diaChADC[clusterChannel2]<{s})))'.format(s=self.saturated_ADC))
			self.trans_grid.canvas['ph{c}_not_saturated'.format(c=c)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
		self.trans_grid.DrawPHGoodAreas('ph{c}_saturated'.format(c='N'), 'clusterCharge{c}'.format(c='N'), '((transparentEvent)&&((diaChADC[clusterChannel0]=={s})||(diaChADC[clusterChannel1]=={s})||(diaChADC[clusterChannel2]=={s})))'.format(s=self.saturated_ADC))
		self.trans_grid.canvas['ph{c}_saturated'.format(c='N')].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift
		self.trans_grid.DrawPHGoodAreas('ph{c}_not_saturated'.format(c='N'), 'clusterCharge{c}'.format(c='N'), '((transparentEvent)&&((diaChADC[clusterChannel0]<{s})&&(diaChADC[clusterChannel1]<{s})&&(diaChADC[clusterChannel2]<{s})))'.format(s=self.saturated_ADC))
		self.trans_grid.canvas['ph{c}_not_saturated'.format(c='N')].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

	def PlotTestForNegative(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_clusters = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr_neg, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr_neg)]
		if cells == 'good':
			list_cuts_noise.append(self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond)
		elif cells == 'bad':
			list_cuts_noise.append(self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
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
			self.trans_grid.DrawProfile2DDiamond('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diaChSignal[clusterChannel{i}]'.format(i=c), '(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]<=-{c})'.format(i=c, c=self.trans_grid.neg_cut))
			# self.trans_grid.DrawProfile2DDiamond('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diaChSignal', '(diaChannels==clusterChannel{i})'.format(i=c))
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
			if len(self.trans_grid.gridAreas.goodAreas_diamond_centers) < 900:
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
				if len(t.trans_grid.gridAreas.goodAreas_diamond_centers) < 900:
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

	def DoAutomatic(self, cells='good', do_save=True):
		self.PlotTestClusterStudies(cells)
		self.PlotTestForNegative(cells)
		self.PlotTest()
		if self.do_saturation:
			self.PlotSaturation()
		if do_save:
			self.SaveCanvas()

	def SetTransparentGrid(self):
		self.trans_grid.SetLines()
		self.trans_grid.CreateTCutGs()
		if self.trans_grid.saturated_ADC != 0:
			self.saturated_ADC = self.trans_grid.saturated_ADC
		else:
			self.saturated_ADC = Get_From_User_Value('saturation_ADC for run ' + str(self.run), 'int', self.trans_grid.saturated_ADC, True)
			self.trans_grid.saturated_ADC = self.saturated_ADC
			self.trans_grid.SavePickle()
		if self.trans_grid.bias != 0:
			self.bias = self.trans_grid.bias
		else:
			self.bias = Get_From_User_Value('bias for run ' + str(self.run), 'float', self.trans_grid.bias, update=True)
			self.trans_grid.bias = self.bias
			self.trans_grid.SavePickle()

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', type='int', default=0, help='run number to be analysed (e.g. 25209)')
	parser.add_option('-a', '--auto', dest='auto', default=False, action='store_true', help='Sets up test, creates plots and saves them automatically if toggled')
	parser.add_option('-c', '--config', dest='config', default='', type='string', help='gives the path to a config file for the test area')

	(options, args) = parser.parse_args()
	run = int(options.run)
	autom = bool(options.auto)
	config = str(options.config)

	t = TestAreas(config, run)
	t.SetTransparentGrid()
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
		t.DoAutomatic('good')

