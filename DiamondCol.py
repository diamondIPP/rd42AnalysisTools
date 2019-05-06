#!/usr/bin/env python
import numpy as np
import ROOT as ro
from Utils import *
from RectangularCell3D import RectangularCell3D
from HexagonalCell3D import HexagonalCell3D

class DiamondCol:
	def __init__(self, coln, numrows, xcenter, lowest_y, cell_height, sides, width_pitch_ratio, run):
		self.col_num = coln
		self.num_rows = numrows
		self.sides = sides
		self.width_pitch_ratio = 1 if self.sides == 4 else width_pitch_ratio
		self.xcenter = xcenter
		self.lowest_y = lowest_y
		self.cell_height = cell_height
		self.run = run
		self.cells = []

	def SetCellsInColumn(self):
		if self.sides == 4:
			self.cells = [RectangularCell3D(self.col_num, ri, self.cell_height, self.run) for ri in xrange(self.num_rows)]
		elif self.sides == 6:
			self.cells = [HexagonalCell3D(self.col_num, ri, self.cell_height, self.width_pitch_ratio, self.run) for ri in xrange(self.num_rows)]
		else:
			ExitMessage('The number of sides is not 4 or 6 but {s}. Cannot create the cells in the column. Exiting'.format(s=self.sides), os.EX_SOFTWARE)
		for ri, cell in enumerate(self.cells):
			cell.SetCellCenter(self.xcenter, self.lowest_y + (2 * ri + 1) * self.cell_height / 2.0)
			cell.CreateTCutG()
			cell.CreateTCutGCenter()

	def GetYCenters(self):
		return np.arange(self.lowest_y + self.cell_height / 2.0, self.lowest_y + (2 * self.num_rows + 1) * self.cell_height / 2.0, self.cell_height, 'f8')

	def GetXYCenters(self):
		ycenters = self.GetYCenters()
		return np.append([ycenters.size * [self.xcenter]], [ycenters], axis=0).T

if __name__ == '__main__':
	z = DiamondCol()