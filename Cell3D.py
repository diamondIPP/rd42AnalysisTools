#!/usr/bin/env python
import numpy as np
import ROOT as ro

class Cell3D:
	def __init__(self, col_num=0, row_num=0, sides=0, run=0):
		self.col_num = col_num
		self.row_num = row_num
		self.sides = sides
		self.run = run
		self.xcenter = 0
		self.ycenter = 0

	def SetCellCenter(self, x, y):
		self.xcenter = x
		self.ycenter = y

	def GetDistanceToCenter(self, x, y):
		return np.sqrt(np.add(np.power(np.subtract(x, self.xcenter, dtype='f8'), 2, dtype='f8'), np.power(np.subtract(y, self.ycenter, dtype='f8'), 2, dtype='f8')), dtype='f8')

if __name__ == '__main__':
	z = Cell3D()