#!/usr/bin/env python

import ROOT as ro
import numpy as np
import ipdb

class GridAreas:
	def __init__(self, numcols=0, numrows=0):
		self.num_cols = numcols
		self.num_rows = numrows
		self.goodAreas_diamond = []
		self.goodAreas_diamond_centers = []
		self.goodAreas_index = []
		self.badAreas_diamond = []
		self.badAreas_diamond_centers = []
		self.badAreas_index = []
		self.goodAreasCutNames_diamond = ''
		self.goodAreasCutNames_diamond_centers = ''
		self.badAreasCutNames_diamond = ''
		self.badAreasCutNames_diamond_centers = ''

	def AddGoodAreas(self, col, row, tcutgs_diamond, tcutgs_diamond_center):
		if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
			self.goodAreas_index.append((col, row))
			tcutgs_diamond[col][row].SetLineColor(ro.kRed)
			tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
			self.goodAreas_diamond.append(tcutgs_diamond[col][row])
			self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddBadAreas(self, col, row, tcutgs_diamond, tcutgs_diamond_center):
		if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
			self.badAreas_index.append((col, row))
			tcutgs_diamond[col][row].SetLineColor(ro.kBlue)
			tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
			self.badAreas_diamond.append(tcutgs_diamond[col][row])
			self.badAreas_diamond_centers.append(tcutgs_diamond_center[col][row])
			tempbad = [cut.GetName() for cut in self.badAreas_diamond]
			self.badAreasCutNames_diamond = '((' + ')||('.join(tempbad) + '))'
			tempbad = [cut.GetName() for cut in self.badAreas_diamond_centers]
			self.badAreasCutNames_diamond_centers = '((' + ')||('.join(tempbad) + '))'

	def AddGoodAreasRow(self, row, coli=0, colf=0, tcutgs_diamond=None, tcutgs_diamond_center=None):
		(colii, colff) = (0, self.num_cols) if coli == 0 and colf == 0 else (coli, colf)
		if 0 <= colii <= colff < self.num_cols and 0 <= row < self.num_rows:
			for col in xrange(colii, colff + 1):
				self.goodAreas_index.append((col, row))
				tcutgs_diamond[col][row].SetLineColor(ro.kRed)
				tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
				self.goodAreas_diamond.append(tcutgs_diamond[col][row])
				self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])

			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddGoodAreasCol(self, col, rowi=0, rowf=0, tcutgs_diamond=None, tcutgs_diamond_center=None):
		(rowii, rowff) = (0, self.num_rows) if rowi == 0 and rowf == 0 else (rowi, rowf)
		if 0 <= col < self.num_cols and 0 <= rowii <= rowff < self.num_rows:
			for row in xrange(rowii, rowff + 1):
				self.goodAreas_index.append((col, row))
				tcutgs_diamond[col][row].SetLineColor(ro.kRed)
				tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
				self.goodAreas_diamond.append(tcutgs_diamond[col][row])
				self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])

			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddRemainingToBadAreas(self, tcutgs_diamond=None, tcutgs_diamond_center=None):
		for col in xrange(0, self.num_cols):
			for row in xrange(0, self.num_rows):
				if (col, row) not in self.goodAreas_index and (col, row) not in self.badAreas_index:
					self.AddBadAreas(col, row, tcutgs_diamond, tcutgs_diamond_center)

	def RemoveFromGoodArea(self, col, row, tcutgs_diamond=None, tcutgs_diamond_center=None):
		if (col, row) in self.goodAreas_index:
			index_g = self.goodAreas_index.index((col, row))
			self.goodAreas_diamond.pop(index_g)
			self.goodAreas_diamond_centers.pop(index_g)
			self.goodAreas_index.pop(index_g)
			self.AddBadAreas(col, row, tcutgs_diamond, tcutgs_diamond_center)
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def ResetAreas(self):
		self.goodAreas_diamond = []
		self.goodAreas_diamond_centers = []
		self.badAreas_diamond = []
		self.badAreas_diamond_centers = []
		self.goodAreas_index = []
		self.badAreas_index = []
		self.goodAreasCutNames_diamond = ''
		self.goodAreasCutNames_diamond_centers = ''
		self.badAreasCutNames_diamond = ''
		self.badAreasCutNames_diamond_centers = ''

if __name__ == '__main__':
	ga = GridAreas(0, 0)
