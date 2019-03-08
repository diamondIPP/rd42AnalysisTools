import numpy as np

class CutManager:
	def __init__(self):
		self.transp_ev = '(transparentEvent)'
		self.tcutg_cell = {}
		self.selected_cells = ''
		self.not_selected_cells = ''
		self.all_cells = '(({s})||({ns}))'
		self.selected_cells_centers = ''
		self.not_selected_cells_centers = ''
		self.no_up_down_borders = '(({l}<=diaChYPred)&&(diaChYPred<={h}))'