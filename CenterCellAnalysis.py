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

		self.suffix = GetSuffixDictionary(self.trans_grid)

		self.phN_chs_var = self.trans_grid.GetPHNChsVar
		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

		self.sat_evts_region_cut = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region_cut = self.trans_grid.cuts_man.not_sat_evts_region

		self.neg_ch_cut = self.trans_grid.cuts_man.GetNegPHChCut
		self.negN_chs_cut = self.trans_grid.cuts_man.GetNegPHNChsCut
		self.not_neg_ch_cut = self.trans_grid.cuts_man.GetNotNegPHChCut
		self.not_negN_chs_cut = self.trans_grid.cuts_man.GetNotNegPHNChsCut

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.in_transp_cluster = ''

		self.sat_evts_region = ''
		self.not_sat_evts_region = ''

		self.noise_varz = {}
		self.noise_friend_varz = {}

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

	def DoPHPlots(self, p_array, cells='all', drawHistos=True, suffixop='', cut='', typ='adc', isFriend=False):
		suffix = self.suffix[cells] if suffixop == '' else suffixop
		for p in p_array:
			if 0 < p < 100:
				varzkey = 'PH{n}_H'.format(n=self.num_strips)
				varz = self.phN_chs_var(self.num_strips, 'H', typ == 'snr', isFriend)
				if p in self.in_central_rect_cut.keys() and p in self.out_central_rect_cut.keys():
					tempcut = cut if cut != '' else self.ph_cuts(varzkey)
					tempcutin = self.trans_grid.cuts_man.AndCuts([tempcut, self.in_central_rect_cut[p]])
					tempcutout = self.trans_grid.cuts_man.AndCuts([tempcut, self.out_central_rect_cut[p]])
					if drawHistos:
						hname = 'PH{n}_H_{t}_cells_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ) if not isFriend else 'PH{n}_H_buffer_{v}_{t}_cells_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, v=self.trans_grid.noise_friend_buffer)
						self.DoCellOverlayPlot(hname, varz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempcutin, cells, typ=typ)
						self.trans_grid.DrawCentralArea(hname, p)
					hname = 'PH{n}_H_{t}_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ) if not isFriend else 'PH{n}_H_buffer_{v}_{t}_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, v=self.trans_grid.noise_friend_buffer)
					self.DoPHPlot(hname, varz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempcutin, cells, typ, drawHistos)

					if drawHistos:
						hname = 'PH{n}_H_{t}_cells_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ) if not isFriend else 'PH{n}_H_buffer_{v}_{t}_cells_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, v=self.trans_grid.noise_friend_buffer)
						self.DoCellOverlayPlot(hname, varz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempcutout, cells, typ=typ)
						self.trans_grid.DrawCentralArea(hname, p)
					hname = 'PH{n}_H_{t}_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ) if not isFriend else 'PH{n}_H_buffer_{v}_{t}_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, v=self.trans_grid.noise_friend_buffer)
					self.DoPHPlot(hname, varz, 'PH{n} highest chs [{t}]'.format(n=self.num_strips, t=typ.upper()), tempcutout, cells, typ, drawHistos)

	def DoCellOverlayPlot(self, name, varz, varname, cut, cells='all', typ='adc'):
		self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname, typ=typ)
		self.PosCanvas(name)

	def DoPHPlot(self, name, varz, varname, cut, cells='all', typ='adc', drawHisto=True):
		self.trans_grid.DrawPHInArea(name, varz, cells, cut, varname=varname, typ=typ, drawHisto=drawHisto)
		if drawHisto: self.PosCanvas(name)

	def DoCenterRegionHistos(self, p_array, cells='all', drawHistos=True, suffixop='', cuts='', typ='adc', isFriend=False):
		self.SetCenterRegionStudies()
		self.DoPHPlots(p_array, cells, drawHistos, suffixop, cuts, typ, isFriend=isFriend)

	def DoCenterRegionStudies(self, dists=np.arange(0.1, 1, 0.025), cells='all', drawHistos=True, suffixop='', cuts='', typ='adc', do_sat=True, isFriend=False):
		self.SetCenterRegionStudies()
		percents = np.unique(np.floor(np.power(dists, 2) * 100 + 0.5).astype('int32'))
		for percent in percents:
			self.trans_grid.CreateTCutGSymmetricRectangle(percent)
		self.DoCenterRegionHistos(percents, cells, drawHistos, suffixop, cuts, typ, isFriend=isFriend)
		self.DrawPHInAndOutGraphs(percents, suffixop, typ, isFriend=isFriend)
		if do_sat: self.DoCenterSaturationStudies(percents, cells, suffixop, cuts, typ)

	def DrawPHInAndOutGraphs(self, percents, suffix, typ='adc', isFriend=False):
		xarray = percents.astype('float64')
		temp0 = np.array([percents[1] - percents[0]] + [percents[i] - percents[i-1] for i in xrange(1, percents.size)])
		temp0fact = (percents[-1] - percents[0]) / (2 * temp0.sum(dtype='float64') - temp0[0] - temp0[-1])
		xarrayerrs = np.multiply(temp0, temp0fact, dtype='float64')

		yarrayin = np.array([self.trans_grid.histo['PH{n}_H{isf}{t}_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')].GetMean() for p in percents], 'float64')
		yarrayinerrs = np.array([self.trans_grid.histo['PH{n}_H{isf}{t}_in_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')].GetMeanError() for p in percents], 'float64')
		yarrayout = np.array([self.trans_grid.histo['PH{n}_H{isf}{t}_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')].GetMean() for p in percents], 'float64')
		yarrayouterrs = np.array([self.trans_grid.histo['PH{n}_H{isf}{t}_out_rect_{p}_percent_{s}'.format(p=p, s=suffix, n=self.num_strips, t=typ, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')].GetMeanError() for p in percents], 'float64')
		maxy = max(yarrayin.max() + yarrayinerrs.max(), yarrayout.max() + yarrayouterrs.max())

		tgraphe_in = ro.TGraphErrors(int(xarray.size), xarray, yarrayin, xarrayerrs, yarrayinerrs)
		ingraphname = 'PH{n}_H{isf}{t}_in_rect_Vs_percent_area_inside_{s}'.format(n=self.num_strips, t=typ, s=suffix, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')
		tgraphe_in.SetNameTitle('g_' + ingraphname, 'g_' + ingraphname)
		tgraphe_in.GetXaxis().SetTitle('percentage of rectangular area inside [%]')
		tgraphe_in.GetYaxis().SetTitle('<PH{n} highest chs inside> [{t}]'.format(n=self.num_strips, t=typ.upper()))
		tgraphe_in.GetYaxis().SetRangeUser(0, maxy)
		tgraphe_in.SetMarkerStyle(8)
		tgraphe_in.SetMarkerColor(ro.kRed)
		tgraphe_in.SetLineColor(ro.kRed)
		self.trans_grid.graph[ingraphname] = tgraphe_in
		self.trans_grid.canvas[ingraphname] = ro.TCanvas('c_' + ingraphname, 'c_' + ingraphname, 1)
		self.trans_grid.graph[ingraphname].Draw('ALP')
		ro.gPad.Update()
		SetDefault1DCanvasSettings(self.trans_grid.canvas[ingraphname])
		self.PosCanvas(ingraphname)

		tgraphe_out = ro.TGraphErrors(int(xarray.size), xarray, yarrayout, xarrayerrs, yarrayouterrs)
		outgraphname = 'PH{n}_H{isf}{t}_out_rect_Vs_percent_area_inside_{s}'.format(n=self.num_strips, t=typ, s=suffix, isf='_' if not isFriend else '_buffer_' + str(self.trans_grid.noise_friend_buffer) + '_')
		tgraphe_out.SetNameTitle('g_' + outgraphname, 'g_' + outgraphname)
		tgraphe_out.GetXaxis().SetTitle('percentage of rectangular area inside [%]')
		tgraphe_out.GetYaxis().SetTitle('<PH{n} highest chs outside> [{t}]'.format(n=self.num_strips, t=typ.upper()))
		tgraphe_out.GetYaxis().SetRangeUser(0, maxy)
		tgraphe_out.SetMarkerStyle(8)
		tgraphe_out.SetMarkerColor(ro.kBlue)
		tgraphe_out.SetLineColor(ro.kBlue)
		self.trans_grid.graph[outgraphname] = tgraphe_out
		self.trans_grid.canvas[outgraphname] = ro.TCanvas('c_' + outgraphname, 'c_' + outgraphname, 1)
		self.trans_grid.graph[outgraphname].Draw('ALP')
		ro.gPad.Update()
		SetDefault1DCanvasSettings(self.trans_grid.canvas[outgraphname])
		self.PosCanvas(outgraphname)

	def DoCenterSaturationStudies(self, percents, cells='all', suffixop='', cuts='', typ='adc'):
		self.SetCenterRegionStudies()
		tempc = self.trans_grid.cuts_man.AndCuts([self.trans_grid.cuts_man.transp_ev, cuts]) if cuts != '' else self.trans_grid.cuts_man.transp_ev
		xarray = percents.astype('float64')
		temp0 = np.array([percents[1] - percents[0]] + [percents[i] - percents[i-1] for i in xrange(1, percents.size)])
		temp0fact = (percents[-1] - percents[0]) / (2 * temp0.sum(dtype='float64') - temp0[0] - temp0[-1])
		xarrayerrs = np.multiply(temp0, temp0fact, dtype='float64')

		ysatin = [self.trans_grid.trans_tree.Draw('event>>hbla', self.trans_grid.cuts_man.AndCuts([tempc, self.sat_adc_inside_cut[p]]), 'goff') for p in percents]
		ysatout = [self.trans_grid.trans_tree.Draw('event>>hbla', self.trans_grid.cuts_man.AndCuts([tempc, self.sat_adc_outside_cut[p]]), 'goff') for p in percents]
		ynosatin = [self.trans_grid.trans_tree.Draw('event>>hbla', self.trans_grid.cuts_man.AndCuts([tempc, self.nosat_adc_inside_cut[p]]), 'goff') for p in percents]
		ynosatout = [self.trans_grid.trans_tree.Draw('event>>hbla', self.trans_grid.cuts_man.AndCuts([tempc, self.nosat_adc_outside_cut[p]]), 'goff') for p in percents]

		ysat_nosat_ratio_in = np.divide(ysatin, np.add(ynosatin, ysatin), dtype='float64')
		ysat_nosat_ratio_out = np.divide(ysatout, np.add(ysatout, ynosatout), dtype='float64')

		ysat_nosat_ratio_in_errs = np.divide(np.sqrt(np.add(np.power(np.multiply(ysatin, np.sqrt(ynosatin)), 2), np.power(np.multiply(ynosatin, np.sqrt(ysatin)), 2))), np.power(np.add(ysatin, ynosatin), 2), dtype='float64')
		ysat_nosat_ratio_out_errs = np.divide(np.sqrt(np.add(np.power(np.multiply(ysatout, np.sqrt(ynosatout)), 2), np.power(np.multiply(ynosatout, np.sqrt(ysatout)), 2))), np.power(np.add(ysatout, ynosatout), 2), dtype='float64')

		tgraphe_sat_ratio_in = ro.TGraphErrors(int(xarray.size), xarray, ysat_nosat_ratio_in, xarrayerrs, ysat_nosat_ratio_in_errs)
		ingraphrationame = 'Sat_events_ratio_in_rect_Vs_percent_area_in_{s}'.format(s=suffixop)
		tgraphe_sat_ratio_in.SetNameTitle('g_' + ingraphrationame, 'g_' + ingraphrationame)
		tgraphe_sat_ratio_in.GetXaxis().SetTitle('percentage of rectangular area inside')
		tgraphe_sat_ratio_in.GetYaxis().SetTitle('sat_tracks/all_tracks')
		tgraphe_sat_ratio_in.GetYaxis().SetRangeUser(0, 1)
		tgraphe_sat_ratio_in.SetMarkerStyle(8)
		tgraphe_sat_ratio_in.SetMarkerColor(ro.kRed)
		tgraphe_sat_ratio_in.SetLineColor(ro.kRed)
		self.trans_grid.graph[ingraphrationame] = tgraphe_sat_ratio_in
		self.trans_grid.canvas[ingraphrationame] = ro.TCanvas('c_' + ingraphrationame, 'c_' + ingraphrationame, 1)
		self.trans_grid.graph[ingraphrationame].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[ingraphrationame])
		self.PosCanvas(ingraphrationame)

		tgraphe_sat_ratio_out = ro.TGraphErrors(int(xarray.size), xarray, ysat_nosat_ratio_out, xarrayerrs, ysat_nosat_ratio_out_errs)
		outgraphrationame = 'Sat_events_ratio_out_rect_Vs_percent_area_in_{s}'.format(s=suffixop)
		tgraphe_sat_ratio_out.SetNameTitle('g_' + outgraphrationame, 'g_' + outgraphrationame)
		tgraphe_sat_ratio_out.GetXaxis().SetTitle('percentage of rectangular area inside')
		tgraphe_sat_ratio_out.GetYaxis().SetTitle('sat_tracks/all_tracks')
		tgraphe_sat_ratio_out.GetYaxis().SetRangeUser(0, 1)
		tgraphe_sat_ratio_out.SetMarkerStyle(8)
		tgraphe_sat_ratio_out.SetMarkerColor(ro.kBlue)
		tgraphe_sat_ratio_out.SetLineColor(ro.kBlue)
		self.trans_grid.graph[outgraphrationame] = tgraphe_sat_ratio_out
		self.trans_grid.canvas[outgraphrationame] = ro.TCanvas('c_' + outgraphrationame, 'c_' + outgraphrationame, 1)
		self.trans_grid.graph[outgraphrationame].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[outgraphrationame])
		self.PosCanvas(outgraphrationame)

	def SetCenterRegionStudies(self):
		self.GetCutsFromCutManager()

	def GetCutsFromCutManager(self):
		self.in_central_rect_cut = self.trans_grid.cuts_man.in_central_rect_region
		self.out_central_rect_cut = self.trans_grid.cuts_man.out_central_rect_region
		self.SetCutsForSaturationRatio()

	def SetCutsForSaturationRatio(self):
		for p, cut in self.in_central_rect_cut.iteritems():
			self.sat_adc_inside_cut[p] = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, cut])
			self.nosat_adc_inside_cut[p] = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, cut])
		for p, cut in self.out_central_rect_cut.iteritems():
			self.sat_adc_outside_cut[p] = self.trans_grid.cuts_man.AndCuts([self.sat_evts_region_cut, cut])
			self.nosat_adc_outside_cut[p] = self.trans_grid.cuts_man.AndCuts([self.not_sat_evts_region_cut, cut])


if __name__ == '__main__':
	c = CenterCellAnalysis(None, 0, 0)
