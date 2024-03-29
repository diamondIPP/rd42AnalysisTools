#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
import cPickle as pickle
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class CellsAnalysis:
	def __init__(self, trans_grid=TransparentGrid()):
		self.window_shift = 3
		self.w = 0
		self.trans_grid = trans_grid
		self.dia_cols = self.trans_grid.dia_cols
		self.numcells = int(self.dia_cols.num_cols * max(self.trans_grid.row_cell_info_diamond['num_even'], self.trans_grid.row_cell_info_diamond['num_odd']))
		self.mean_ph_cell_dic = {typ: {} for typ in ['adc', 'snr']}
		self.mean_ph_cell_friend_dic = {typ: {} for typ in ['adc', 'snr']}
		self.suffix = GetSuffixDictionary(self.trans_grid)
		self.qvalues = {cell: {} for cell in ['', 'all', 'bad', 'good']}

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def GetMeanPHWithinCell(self, col, row, var='clusterChargeN', bins=720, liminf=-2400, limsup=4800, cut='', corner=-1):
		"""
		Method that returns the mean PH of a variable 'var' of type 'typ' within the cell in column 'col' and row 'row'
		:param var: variable to calculate the mean
		:param bins: bins to use in the histogram to calculate the mean
		:param liminf: histogram's inferior limit to calculate the mean
		:param limsup: histogram's superior limit to calculate the mean
		:param cut: extra cuts besides the spatial location within the cell
		:param corner: if corner is 0, then the region within the read-out column is taken into account. If 0 < corner <= num_sides,
						then the region surrounding one of the bias columns within a cell is taken into account. If corner is another value,
						then the whole cell is taken into account
		:return: returns the mean value of the histogram calculated with the specified parameters
		"""
		if corner == 0:
			tempc = self.trans_grid.cuts_man.AndCuts(['transparentEvent', self.dia_cols.cols[col].cells[row].cutg_read_out, cut])
		elif 0 < corner <= self.trans_grid.num_sides:
			tempc = self.trans_grid.cuts_man.AndCuts(['transparentEvent', self.dia_cols.cols[col].cells[row].cutg_bias[corner - 1], cut])
		else:
			tempc = self.trans_grid.cuts_man.AndCuts(['transparentEvent', self.dia_cols.cols[col].cells[row].cutg.GetName(), cut])
		draw_string = '{v}>>temphrc({nb},{lb},{hb})'.format(v=var, nb=bins, lb=liminf, hb=limsup)
		self.trans_grid.trans_tree.Draw(draw_string, tempc, 'goff')
		temph = ro.gDirectory.Get('temphrc')
		return temph.GetMean()

	def LoadMeanPHPerCell(self, typ='adc', isFriend=False):
		"""
		Loads the whole dictionary with the information stored in the pickle file of type 'typ'
		:param typ: specifies whether to load the pickle with 'adc' information, or the one with 'snr'
		:param isFriend: if True, it loads the information coming from the pickle that has the data calculated from a pedTree friend
		:return:
		"""
		pathp = '{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir, t=typ) if not isFriend else '{d}/{r}/{s}/mean_ph_cell_friend_{t}.{r}.pkl'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir, t=typ)
		if os.path.isfile(pathp):
			with open(pathp) as pkl:
				means_temp = pickle.load(pkl)
				if means_temp:
					if not isFriend:
						self.mean_ph_cell_dic[typ] = means_temp
					else:
						self.mean_ph_cell_friend_dic[typ] = means_temp

	def IsVariableInMeanPHDic(self, var='clusterChargeN', typ='adc', isFriend=False):
		"""
		This method returns True or False if there is info on the mean ph of a cell for a certain type, variable and data set
		:param var: Is a string describing the PH variable requested e.g. clusterChargeN
		:param typ: String which is either 'adc' or 'snr'
		:param isFriend: If true, the information is from a friend pedTree which has different pedestal and therefore ph calculations
		:return: boolean if the variable is found within the corresponding dictionary
		"""
		if not isFriend:
			if var in self.mean_ph_cell_dic[typ].keys():
				return True
		else:
			if var in self.mean_ph_cell_friend_dic[typ].keys():
				return True
		return False

	def CalculateCellMeanPH(self, var='clusterChargeN', typ='adc', isFriend=False):
		"""
		This method calculates the mean of a variable 'var' of type 'typ' on a dataset that could be from a pedTree friend, and stores it
		in the corresponding dictionary
		:param var: variable used to calculate the meanPH
		:param typ: specifies if the variable is of type 'adc' or 'snr'
		:param isFriend: if true, the results come from a pedTree friend. the variable has to be consistent
		:return:
		"""
		tempbar = CreateProgressBarUtils(self.numcells)
		tempbar.start()
		if not isFriend:
			self.mean_ph_cell_dic[typ][var] = {}
		else:
			self.mean_ph_cell_friend_dic[typ][var] = {}
		for col in xrange(self.dia_cols.num_cols):
			if not isFriend:
				self.mean_ph_cell_dic[typ][var][col] = {}
			else:
				self.mean_ph_cell_friend_dic[typ][var][col] = {}
			for row in xrange(self.dia_cols.cols[col].num_rows):
				temphmean = self.GetMeanPHWithinCell(col, row, var, 720, -2400 if typ == 'adc' else -240, 4800 if typ == 'adc' else 480)
				if not isFriend:
					self.mean_ph_cell_dic[typ][var][col][row] = temphmean
				else:
					self.mean_ph_cell_friend_dic[typ][var][col][row] = temphmean
				tempbar.update(col * self.dia_cols.cols[col].num_rows + row + 1)
		tempbar.finish()

	def SavePHMeanDicPickle(self, typ='adc', isFriend=False):
		"""
		Saves the mean PH for each cell in a pickle file of type 'typ'
		:param typ: indicates if the information to save is of type 'adc' or 'snr'
		:param isFriend: if True, it indicates that the information to save was calculated from a pedTree friend.
		:return:
		"""
		meanphobj = self.mean_ph_cell_dic[typ] if not isFriend else self.mean_ph_cell_friend_dic[typ]
		if not os.path.isdir('{d}/{r}/{s}'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir)):
			os.makedirs('{d}/{r}/{s}'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir))
		print 'Saving mean PH of type', typ, 'in file', '{d}/{r}/{s}'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir), '...', ; sys.stdout.flush()
		if not isFriend:
			pickle.dump(meanphobj, open('{d}/{r}/{s}/mean_ph_cell_{t}.{r}.pkl'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir, t=typ), 'wb'))
		else:
			pickle.dump(meanphobj, open('{d}/{r}/{s}/mean_ph_cell_friend_{t}.{r}.pkl'.format(d=self.trans_grid.dir, r=self.trans_grid.run, s=self.trans_grid.pkl_sbdir, t=typ), 'wb'))
		print 'Done'

	def FindCellsQuality(self, cells='all'):
		print 'Get relative charge ratio of neighbor channels...'
		hname = self.Draw2DChargeRatioWithNeighborStrips(cells, suffix='{s}_no_cuts'.format(s=self.suffix[cells]))
		# Reduced the bin content to ignore the bins with lower than 5% of the histogram Maximum
		threshold = self.trans_grid.histo[hname].GetMaximum() / 20.
		histor = ReduceHistoContent(self.trans_grid.histo[hname], threshold, True)
		hname_reduced = histor.GetName()[2:]
		self.trans_grid.histo[hname_reduced] = histor
		print 'Project in vertical axis...'
		# Project reduced histogram in vertical axis to find the peak(s) of the distribution
		historpy = histor.ProjectionY('h_' + hname_reduced + '_py', 0, -1, 'e')
		self.trans_grid.histo[hname_reduced + '_py'] = historpy
		historpy.SetTitle('h_' + hname_reduced + '_py')
		historpy.GetXaxis().SetRangeUser(historpy.GetMean() - 2 * historpy.GetRMS(), historpy.GetMean() + 2 * historpy.GetRMS())
		self.trans_grid.canvas[hname_reduced + '_py'] = ro.TCanvas('c_' + hname_reduced + '_py', 'c_' + hname_reduced + '_py', 1)
		historpy.Draw('e hist')
		SetDefault1DStats(historpy)
		self.PosCanvas(hname_reduced + '_py')
		self.trans_grid.FitTwoGaus(hname_reduced + '_py', 2)
		fitpy = self.trans_grid.fits[hname_reduced + '_py']
		fitpyparams = fitpy.GetParams()
		fitpyparams = np.array([fitpyparams[it] for it in xrange(6)], 'f8')
		if fitpyparams[1] < fitpyparams[4]:
			muy1 = fitpyparams[1]
			sy1 = fitpyparams[2]
			muy2 = fitpyparams[4]
			sy2 = fitpyparams[5]
		else:
			muy1 = fitpyparams[4]
			sy1 = fitpyparams[5]
			muy2 = fitpyparams[1]
			sy2 = fitpyparams[2]
		histor.GetYaxis().SetRangeUser(muy1 - 2. * sy1, muy2 + 2. * sy2)
		print 'Project in horizontal axis...'
		# Project reduced histogram around the lowest gaussian peak on the horizontal axis
		historpx = histor.ProjectionX('h_' + hname_reduced + '_px', 0, -1, 'e')
		self.trans_grid.histo[hname_reduced + '_px'] = historpx
		historpx.SetTitle('h_' + hname_reduced + '_px')
		historpx.GetXaxis().SetRangeUser(historpx.GetMean() - 2 * historpx.GetRMS(), historpx.GetMean() + 2 * historpx.GetRMS())
		self.trans_grid.canvas[hname_reduced + '_px'] = ro.TCanvas('c_' + hname_reduced + '_px', 'c_' + hname_reduced + '_px', 1)
		historpx.Draw('e hist')
		self.PosCanvas(hname_reduced + '_px')
		SetDefault1DStats(historpx)
		self.trans_grid.FitTwoGaus(hname_reduced + '_px', 2)
		fitpx = self.trans_grid.fits[hname_reduced + '_px']
		fitpxparams = fitpx.GetParams()
		fitpxparams = np.array([fitpxparams[it] for it in xrange(6)], 'f8')
		if fitpxparams[2] < fitpxparams[5]:
			mux1 = fitpxparams[1]
			sx1 = fitpxparams[2]
			mux2 = fitpxparams[4]
			sx2 = fitpxparams[5]
		else:
			mux1 = fitpxparams[4]
			sx1 = fitpxparams[5]
			mux2 = fitpxparams[1]
			sx2 = fitpxparams[2]
		print 'Setting cuts...', ; sys.stdout.flush()
		# Create Cuts to select good and bad regions in the ph ch ratio space
		self.trans_grid.cuts_man.SetQualityRatioCuts(mux1, sx1, mux2, sx2, muy1, sy1, muy2, sy2, varch0=self.trans_grid.GetPHChVar(0, 'Ch'), varch1=self.trans_grid.GetPHChVar(1, 'Ch'), varch2=self.trans_grid.GetPHChVar(2, 'Ch'))
		print 'Done'
		print 'Plot events according to cuts...', ; sys.stdout.flush()
		hname = self.Draw2DChargeRatioWithNeighborStrips(cells, cuts=self.trans_grid.cuts_man.cut_good_ratio_region, suffix='{s}_good_region'.format(s=self.suffix[cells]))
		self.PosCanvas(hname)
		hname = self.Draw2DChargeRatioWithNeighborStrips(cells, cuts=self.trans_grid.cuts_man.cut_bad_ratio_region, suffix='{s}_bad_region'.format(s=self.suffix[cells]))
		self.PosCanvas(hname)
		print 'Done'
		print 'Create charge maps of tracks inside or outside the created regions in the channel pulse height ratio space...'
		# Create charge maps of tracks that are inside the good and bad regions in the ph ch ratio space
		self.DrawRatioTracksFromQualityCuts(cells, suffix='{s}_no_cuts'.format(s=self.suffix[cells]))

	def Draw2DChargeRatioWithNeighborStrips(self, cells='all', cuts='', suffix=''):
		deltax, deltay = 0.004, 0.004
		hlimsx, hlimsy = Get1DLimits(-1, 5, deltax), Get1DLimits(-1, 3, deltay)
		hname = 'relative_charge_ratio_neighbor_chs'
		hname = hname if suffix == '' else hname + '_' + suffix
		tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cuts, cells)
		self.trans_grid.DrawHisto2D(hname, hlimsx['min'], hlimsx['max'], deltax, 'PH_Ch1/PH_Ch0', hlimsy['min'], hlimsy['max'], deltay, 'PH_Ch2/PH_Ch0', '{n1}/{n0}'.format(n1=self.trans_grid.GetPHChVar(1, 'Ch'), n0=self.trans_grid.GetPHChVar(0, 'Ch')), '{n2}/{n0}'.format(n2=self.trans_grid.GetPHChVar(2, 'Ch'), n0=self.trans_grid.GetPHChVar(0, 'Ch')), tempc)
		self.trans_grid.canvas[hname].SetGridx()
		self.trans_grid.canvas[hname].SetGridy()
		self.PosCanvas(hname)
		return hname

	def DrawRatioTracksFromQualityCuts(self, cells='all', cuts='', suffix=''):
		self.DrawChargeMapsTracksFromQualityCuts(cells, cuts, suffix)
		hname_good = 'PH_Ch0_tracks_good_ph_ch_ratio_region'
		hname_good += '' if suffix == '' else '_' + suffix
		hname_bad = 'PH_Ch0_tracks_bad_ph_ch_ratio_region'
		hname_bad += '' if suffix == '' else '_' + suffix
		# Calculate ratio between well behaved tracks and bad behaved tracks
		print 'Calculate ratio of number of tracks inside and outside the cuts'
		hname_ratio = 'ratio_tracks_good_over_bad_ph_ch_ratio_region'
		hname_ratio += '' if suffix == '' else '_' + suffix
		hname_denom = 'tracks_bad_ph_ch_ratio_region'
		hname_denom += '' if suffix == '' else '_' + suffix
		histo_ratio = self.trans_grid.histo['hit_map_' + hname_good].Clone('h_' + hname_ratio)
		histo_ratio.SetTitle('h_' + hname_ratio)
		histo_denom = self.trans_grid.histo['hit_map_' + hname_bad].Clone('h_' + hname_denom)
		histo_denom.SetTitle('h_' + hname_denom)
		# for taking into account empty cells:
		for bini in xrange(1, histo_ratio.GetNbinsX() * histo_ratio.GetNbinsY() + 1):
			histo_ratio.AddBinContent(bini, 1e-6)
		for bini in xrange(1, histo_denom.GetNbinsX() * histo_denom.GetNbinsY() + 1):
			histo_denom.AddBinContent(bini, 1e-6)
		histo_ratio.Divide(histo_denom)
		histo_ratio.SetMinimum(0)
		histo_ratio.SetMaximum(2)
		for bini in xrange(1, histo_ratio.GetNbinsX() * histo_ratio.GetNbinsY() + 1):
			binic = histo_ratio.GetBinContent(bini)
			if binic <= 2:
				continue
			else:
				histo_ratio.SetBinContent(bini, 2)
		# Show results...
		print 'Plotting results...', ; sys.stdout.flush()
		self.trans_grid.canvas[hname_ratio] = ro.TCanvas('c_' + hname_ratio, 'c_' + hname_ratio, 1)
		self.trans_grid.histo[hname_ratio] = histo_ratio
		histo_ratio.Draw('colz')
		# self.trans_grid.colorPalleteExec1.Draw()
		# histo_ratio.Draw('colz same')
		self.trans_grid.DrawTCutGs(hname_ratio, 'diamond')
		self.PosCanvas(hname_ratio)
		hname2 = 'cells_quality'
		hname2 += '' if suffix == '' else '_' + suffix
		histo_q2 = histo_ratio.Clone('h_' + hname2)
		histo_q2.SetTitle('h_' + hname2)
		histo_q2.Reset()
		histo_q2.GetZaxis().SetTitle('quality [a.u.]')
		self.qvalues[cells] = {}
		for col in self.trans_grid.dia_cols.cols:
			self.qvalues[cells][col.col_num] = {}
			for cell in col.cells:
				self.qvalues[cells][col.col_num][cell.row_num] = GetAverageInCell(histo_ratio, cell.cutg)
		bx, by, bz = np.zeros(1, 'i4'), np.zeros(1, 'i4'), np.zeros(1, 'i4')
		x, y, z = 0, 0, 0
		for bini in xrange(1, histo_q2.GetNbinsX() * histo_q2.GetNbinsY() + 1):
			histo_q2.GetBinXYZ(bini, bx, by, bz)
			x = histo_q2.GetXaxis().GetBinCenter(int(bx))
			y = histo_q2.GetYaxis().GetBinCenter(int(by))
			[col, row] = self.dia_cols.GetColNumRowNumOfPoint(x, y)
			if [col, row] != [-999, -999]:
				histo_q2.SetBinContent(bini, self.qvalues[cells][col][row])
		self.trans_grid.canvas[hname2] = ro.TCanvas('c_' + hname2, 'c_' + hname2, 1)
		self.trans_grid.histo[hname2] = histo_q2
		histo_q2.Draw('colz')
		# self.trans_grid.colorPalleteExec1.Draw()
		# histo_q2.Draw('colz same')
		SetDefault2DStats(histo_q2)
		self.trans_grid.DrawTCutGs(hname2, 'diamond')
		self.PosCanvas(hname2)
		histo_q2.SetMaximum()
		ro.gPad.RedrawAxis()
		# ro.gStyle.SetPalette(55)
		# ro.gStyle.SetNumberContours(50)
		hname1 = 'cells_quality_distribution'
		hname1 += '' if suffix == '' else '_' + suffix
		ncells = self.dia_cols.num_cols * max(self.trans_grid.row_cell_info_diamond['num_even'], self.trans_grid.row_cell_info_diamond['num_odd'])
		histo1 = ro.TH1F('h_' + hname1, 'h_' + hname1, ncells, 0, 2)
		self.trans_grid.histo[hname1] = histo1
		for value in [val for qcv in self.qvalues[cells].values() for val in qcv.values()]:
			histo1.Fill(value)
		self.trans_grid.canvas[hname1] = ro.TCanvas('c_' + hname1, 'c_' + hname1, 1)
		histo1.Draw('e hist')
		SetDefault1DStats(histo1)
		SetDefault1DCanvasSettings(self.trans_grid.canvas[hname1])
		self.PosCanvas(hname1)
		histo1c = histo1.GetCumulative()
		histo1c.Sumw2()
		histo1c.Scale(1.0 / histo1c.GetMaximum())
		self.trans_grid.histo[hname1 + '_cumulative'] = histo1c
		self.trans_grid.canvas[hname1 + '_cumulative'] = ro.TCanvas('c_{n}_cumulative'.format(n=hname1), 'c_{n}_cumulative'.format(n=hname1), 1)
		histo1c.Draw('hist')
		SetDefault1DStats(histo1c)
		SetDefault1DCanvasSettings(self.trans_grid.canvas[hname1 + '_cumulative'])
		self.PosCanvas(hname1 + '_cumulative')
		print 'Done'

	def DrawChargeMapsTracksFromQualityCuts(self, cells='all', cuts='', suffix=''):
		if self.trans_grid.cuts_man.cut_good_ratio_region != '' and self.trans_grid.cuts_man.cut_bad_ratio_region != '':
			hname_good = 'PH_Ch0_tracks_good_ph_ch_ratio_region'
			hname_good += '' if suffix == '' else '_' + suffix
			hname_bad = 'PH_Ch0_tracks_bad_ph_ch_ratio_region'
			hname_bad += '' if suffix == '' else '_' + suffix
			tempcg = self.trans_grid.cuts_man.AndCuts([self.trans_grid.cuts_man.cut_good_ratio_region, cuts])
			tempcb = self.trans_grid.cuts_man.AndCuts([self.trans_grid.cuts_man.cut_bad_ratio_region, cuts])
			self.trans_grid.DrawProfile2DDiamond(hname_good, self.trans_grid.GetPHChVar(0, 'Ch'), 'PH cluster channel 0 [ADC]', cells, tempcg)
			self.PosCanvas(hname_good)
			self.trans_grid.DrawProfile2DDiamond(hname_bad, self.trans_grid.GetPHChVar(0, 'Ch'), 'PH cluster channel 0 [ADC]', cells, tempcb)
			self.PosCanvas(hname_bad)
			self.trans_grid.GetOccupancyFromProfile(hname_good)
			self.PosCanvas('hit_map_' + hname_good)
			self.trans_grid.GetOccupancyFromProfile(hname_bad)
			self.PosCanvas('hit_map_' + hname_bad)
		else:
			print 'Can\'t do the plots. Must first calculate cut_good_ratio_region and cut_bad_ratio_region which are generated by cut manager. Run First FindCellsQuality in transparent grid'

	def Do3DColumnsStudies(self, cells='all'):
		suffix = '{s}_all_cell'.format(s=self.suffix[cells])
		self.Draw2DChargeRatioWithNeighborStrips(cells, '', suffix)
		suffix = '{s}_RO_columns'.format(s=self.suffix[cells])
		self.Draw2DChargeRatioWithNeighborStrips(cells, self.trans_grid.cuts_man.ro_3Dcol_cut, suffix)
		suffix = '{s}_bias_columns'.format(s=self.suffix[cells])
		self.Draw2DChargeRatioWithNeighborStrips(cells, self.trans_grid.cuts_man.bias_all_3Dcol_cut, suffix)




if __name__ == '__main__':
	c = CellsAnalysis(TransparentGrid())
