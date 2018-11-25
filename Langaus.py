#!/usr/bin/env python

import ROOT as ro
import numpy as np
import ipdb

class LanGaus:
	def __init__(self, histo):
		self.params = np.zeros(4, 'f8')
		self.conv_steps = 1000
		self.sigma_conv = 5
		self.mpshift = -0.22278298
		self.paramErrors = np.zeros(4, 'f8')
		self.paramsLimitsLow = np.zeros(4, 'f8')
		self.paramsLimitsHigh = np.zeros(4, 'f8')
		self.paramsFitErrors = np.zeros(4, 'f8')
		self.histo = histo
		self.chi2 = 0
		self.ndf = 0
		self.fit_range = [histo.GetBinLowEdge(histo.FindFirstBinAbove()), histo.GetBinLowEdge(histo.FindLastBinAbove()) + histo.GetBinWidth(histo.FindLastBinAbove())]
		self.fit = None
		self.EstimateParameters()
		self.EstimateParametersLimits()

	def LangausFunc(self, x, params):
		mpc = params[1] - self.mpshift * params[0]
		xlow, xhigh = [x[0] + self.sigma_conv * i * params[3] for i in [-1, 1]]
		step = (xhigh - xlow) / self.conv_steps
		sums = 0
		i = 1
		while i <= int(np.ceil(self.conv_steps/2.0)):
			xx = xlow + (i - 0.5) * step
			fland = ro.TMath.Landau(xx, mpc, params[0]) / params[0]
			sums += fland + ro.TMath.Gaus(x[0], xx, params[3])
			xx = xhigh - (i - 0.5) * step
			fland = ro.TMath.Landau(xx, mpc, params[0]) / params[0]
			sums += fland * ro.TMath.Gaus(x[0], xx, params[3])
			i += 1
		return params[2] * step * sums / (np.sqrt(2 * np.pi, dtype='f8') * params[3])


	def EstimateParameters(self):
		self.params[0] = self.histo.GetRMS()/5.74
		self.params[1] = ((self.histo.GetMean() + self.histo.GetBinCenter(self.histo.GetMaximumBin())) / 2.0) + self.mpshift * self.params[0]
		self.params[2] = self.histo.Integral() * 10
		self.params[3] = self.histo.GetRMS() * 3.55

	def EstimateParametersLimits(self):
		self.paramsLimitsLow[0] = self.histo.GetRMS()/10.0
		self.paramsLimitsLow[1] = max(0, self.histo.GetBinCenter(self.histo.GetMaximumBin()) + self.mpshift * 2*self.params[0])
		self.paramsLimitsLow[2] = self.histo.Integral() * 5.0
		self.paramsLimitsLow[3] = self.histo.GetRMS() * 1.5
		self.paramsLimitsHigh[0] = self.histo.GetRMS()/4.0
		self.paramsLimitsHigh[1] = min(self.histo.GetXaxis().GetXmax(), self.histo.GetBinCenter(self.histo.GetMaximumBin()) - self.mpshift * 2*self.params[0])
		self.paramsLimitsHigh[2] = self.histo.Integral() * 15
		self.paramsLimitsHigh[3] = self.histo.GetRMS()*5.0

	def LanGausFit(self, nconv=1000):
		fit_name = 'fit_{n}'.format(n=self.histo.GetName())
		self.conv_steps = nconv

		fit_old = ro.gROOT.GetListOfFunctions().FindObject(fit_name)
		if fit_old:
			fit_old.Delete()
			del fit_old
		self.fit = ro.TF1(fit_name, self.LangausFunc, self.histo.GetXaxis().GetXmin(), self.histo.GetXaxis().GetXmax(), 4)
		self.fit.SetNpx(1000)
		self.fit.SetParameters(self.params)
		self.fit.SetParNames('Width', 'MP', 'Area', 'GSigma')

		for i in xrange(len(self.params)):
			self.fit.SetParLimits(i, self.paramsLimitsLow[i], self.paramsLimitsHigh[i])
		ipdb.set_trace()
		self.histo.Fit(fit_name, 'QRB0', '', self.fit_range[0], self.fit_range[1])
		self.fit.GetParameters(self.params)
		for i in xrange(len(self.params)):
			self.paramsFitErrors[i] = self.fit.GetParError(i)
		self.chi2 = self.fit.GetChisquare()
		self.ndf = self.fit.GetNDF()
