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


if __name__ == '__main__':
	c = CellsAnalysis(TransparentGrid())
