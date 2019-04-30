#!/usr/bin/env python
import numpy as np
import ROOT as ro
from Utils import *
from DiamondCol import DiamondCol

class DiamondColumns:
	def __init__(self, num_cols, cell_height, sides, run):
		self.num_cols = num_cols
		self.cell_height = cell_height
		self.sides = sides
		self.run = run
		self.cols = []
		# self.cols = {}

	def SetupColumns(self, col, numrows, xcenter, lowest_y):
		coli = int(RoundInt(col))
		self.cols.append(DiamondCol(coli, numrows, xcenter, lowest_y, self.cell_height, self.sides, self.run))
		# self.cols[coli] = DiamondCol(coli, numrows, xcenter, lowest_y, self.cell_height, self.sides, self.run)

	def GetColYCenters(self, col):
		return self.cols[col].GetYCenters()

	def GetColXCenter(self, col):
		return self.cols[col].xcenter

if __name__ == '__main__':
	z = DiamondColumns()