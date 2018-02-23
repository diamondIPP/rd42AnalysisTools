#!/usr/bin/env python
from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, kBlue, kBlack, gStyle, TBrowser, Long, TF1, gROOT, TGraph, TMultiGraph, kMagenta
from optparse import OptionParser
from numpy import array, floor
from copy import deepcopy
import os, sys

__author__ = 'DA'

class SecondPulseID:
    def __init__(self, run=22011, sourceDir='./', outputDir='', low=1500, top=1700):
        print 'Creating SecondPulseID instance for run:', run
        self.run = run
        self.sourceDir = sourceDir
        self.outputDir = outputDir if outputDir != '' else '{s}/{r}/pedestalAnalysis/'.format(r=self.run, s=self.sourceDir)
        self.ADClow = low
        self.ADChigh = top
        self.selRootFile = TFile(sourceDir + '/{r}/selectionData.{r}.root'.format(r=self.run))
        self.selTree = self.selRootFile.Get('selectionTree')
        self.vecSize = {'ADC': self.selTree.GetBranch('DiaADC').GetLeaf('DiaADC').GetLen(), 'Ped': self.selTree.GetBranch('diaPedestalMean').GetLeaf('diaPedestalMean').GetLen(), 'CMN': self.selTree.GetBranch('commonModeNoise').GetLeaf('commonModeNoise').GetLen(), 'Signal': self.selTree.GetBranch('DiaADC').GetLeaf('DiaADC').GetLen()}
        self.adc, self.ped, self.cmn, self.event, self.sig, self.ch = [], [], [], [], [], []
        self.adcEvent, self.pedEvent, self.cmnEvent, self.sigEvent, self.chEvent, self.adcCh, self.pedCh, self.cmnCh, self.sigCh = TGraph(), TGraph(), TGraph(), TGraph(), TGraph(), TGraph(), TGraph(), TGraph(), TGraph()
        self.multiEvent = TMultiGraph()
        self.multiCh = TMultiGraph()

        print 'Loading ADC, Pedestal and cmn vectors...', ; sys.stdout.flush()
        self.BranchesToVectors()
        print 'Done'

        print 'Creating TGraphs...', ; sys.stdout.flush()
        self.VectorsToGraphs()
        print 'Done'

        # print 'Loading ADC vector... ', ; sys.stdout.flush()
        # self.DiaADCVectorEventCh = self.BranchToVector('DiaADC')
        # print 'Done'
        # print 'Loading Pedestals vector... ', ; sys.stdout.flush()
        # self.DiaPedVectorEventCh = self.BranchToVector('diaPedestalMean')
        # print 'Done'
        # print 'Loading CMN vector... ', ; sys.stdout.flush()
        # self.DiaCMNVectorEvent = self.BranchToVector('commonModeNoise')
        # print 'Done'
        # print 'Loading cmn Correction flag vector...', ; sys.stdout.flush()
        # self.DiaCMNCorrVectEvent = self.BranchToVector('cmnCorrection')
        # print 'Done'

        # self.noiseRootFile = TFile(sourceDir + '/{r}/pedestalAnalysis/root/hNoiseDistributionOfAllNonHitChannels_DiaChannel.{r}.root'.format(r=self.run))
        # self.noiseCMnootFile = TFile(sourceDir + '/{r}/pedestalAnalysis/root/hNoiseDistributionOfAllNonHitChannels_DiaChannel_CMNcorrected.{r}.root'.format(r=self.run))
        # self.noiseHisto2D = TH2F(self.noiseRootFile.Get('cRoot_hNoiseDistributionOfAllNonHitChannels_DiaChannel').GetPrimitive('hNoiseDistributionOfAllNonHitChannels_DiaChannel'))
        # self.noiseCMNHisto2D = TH2F(self.noiseCMnootFile.Get('cRoot_hNoiseDistributionOfAllNonHitChannels_DiaChannel_CMNcorrected').GetPrimitive('hNoiseDistributionOfAllNonHitChannels_DiaChannel_CMNcorrected'))
        # self.noiseHisto2D.GetXaxis().SetTitle('Diamond ch')
        # self.noiseHisto2D.GetYaxis().SetTitle('adc-ped')
        # self.noiseHisto2D.GetZaxis().SetTitle('# entries')
        # self.noiseHisto2D.GetYaxis().SetRangeUser(-32,32)
        # self.noiseCMNHisto2D.GetXaxis().SetTitle('Diamond ch')
        # self.noiseCMNHisto2D.GetYaxis().SetTitle('adc-ped-cmn')
        # self.noiseCMNHisto2D.GetZaxis().SetTitle('# entries')
        # self.noiseCMNHisto2D.GetYaxis().SetRangeUser(-32,32)
        # self.noiseHisto1D = 0
        # self.noiseCMNHisto1D = 0
        # self.diamondNoiseChannels = {'ch_low': 74, 'ch_high': 80}
        # self.nameChSelection = 'Strip'
        # self.fidcut = 0
        # self.fidcutCMN = 0
        # self.noiseHisto2D_fid = 0
        # self.noiseCMNHisto2D_fid = 0
        # self.fitNoise = 0
        # self.fitNoiseCMNC = 0
        # if not os.path.isdir('{dir}/Plots'.format(dir=self.outputDir)):
        #     os.makedirs('{dir}/Plots'.format(dir=self.outputDir))
        # if not os.path.isdir('{dir}/root'.format(dir=self.outputDir)):
        #     os.makedirs('{dir}/root'.format(dir=self.outputDir))
        # gStyle.SetPalette(55)
        # gStyle.SetNumberContours(999)
        self.bla = []

    def BranchesToVectors(self):
        leng = self.selTree.Draw('DiaADC:diaPedestalMean:commonModeNoise:EventNumber','({l}<DiaADC-diaPedestalMean-commonModeNoise)&&(DiaADC-diaPedestalMean-commonModeNoise<{h})&&(IsInFiducialCut)&&(hasValidSiliconTrack&&fiducialRegion==3)'.format(l=self.ADClow, h=self.ADChigh),'goff')
        if leng > 1000000:
            self.selTree.SetEstimate(leng)
            leng = self.selTree.Draw('DiaADC:diaPedestalMean:commonModeNoise:EventNumber','({l}<DiaADC-diaPedestalMean-commonModeNoise)&&(DiaADC-diaPedestalMean-commonModeNoise<{h})&&(IsInFiducialCut)&&(hasValidSiliconTrack&&fiducialRegion==3)'.format(l=self.ADClow, h=self.ADChigh),'goff')
        self.adc = [self.selTree.GetV1()[ev] for ev in xrange(leng)]
        self.ped = [self.selTree.GetV2()[ev] for ev in xrange(leng)]
        self.cmn = [self.selTree.GetV3()[ev] for ev in xrange(leng)]

        leng = self.selTree.Draw('EventNumber:diaChannel','({l}<DiaADC-diaPedestalMean-commonModeNoise)&&(DiaADC-diaPedestalMean-commonModeNoise<{h})&&(IsInFiducialCut)&&(hasValidSiliconTrack&&fiducialRegion==3)'.format(l=self.ADClow, h=self.ADChigh),'goff')
        if leng > 1000000:
            self.selTree.SetEstimate(leng)
            leng = self.selTree.Draw('EventNumber:diaChannel','({l}<DiaADC-diaPedestalMean-commonModeNoise)&&(DiaADC-diaPedestalMean-commonModeNoise<{h})&&(IsInFiducialCut)&&(hasValidSiliconTrack&&fiducialRegion==3)'.format(l=self.ADClow, h=self.ADChigh),'goff')

        self.event = [self.selTree.GetV1()[ev] for ev in xrange(leng)]
        self.ch = [self.selTree.GetV2()[ev] for ev in xrange(leng)]
        self.sig = [self.adc[ev]-self.ped[ev]-self.cmn[ev] for ev in xrange(leng)]

    def VectorsToGraphs(self):
        self.adcEvent = TGraph(len(self.adc), array(self.event,'d'), array(self.adc,'d'))
        self.adcEvent.SetLineColor(kBlack)
        self.pedEvent = TGraph(len(self.ped), array(self.event,'d'), array(self.ped,'d'))
        self.pedEvent.SetLineColor(kRed)
        self.cmnEvent = TGraph(len(self.cmn), array(self.event,'d'), array(self.cmn,'d'))
        self.cmnEvent.SetLineColor(kMagenta)
        self.sigEvent = TGraph(len(self.sig), array(self.event,'d'), array(self.sig,'d'))
        self.sigEvent.SetLineColor(kBlue)

        self.adcCh = TGraph(len(self.adc), array(self.ch,'d'), array(self.adc,'d'))
        self.adcCh.SetLineColor(kBlack)
        self.pedCh = TGraph(len(self.ped), array(self.ch,'d'), array(self.ped,'d'))
        self.pedCh.SetLineColor(kRed)
        self.cmnCh = TGraph(len(self.cmn), array(self.ch,'d'), array(self.cmn,'d'))
        self.cmnCh.SetLineColor(kMagenta)
        self.sigCh = TGraph(len(self.sig), array(self.ch,'d'), array(self.sig,'d'))
        self.sigCh.SetLineColor(kBlue)

        self.chEvent = TGraph(len(self.ch), array(self.event,'d'), array(self.ch,'d'))
        self.chEvent.SetLineColor(kBlack)

        self.multiEvent.SetNameTitle('MultiGraph_event', 'MultiGraph_event')
        self.multiEvent.Add(self.adcEvent)
        self.multiEvent.Add(self.pedEvent)
        self.multiEvent.Add(self.cmnEvent)
        self.multiEvent.Add(self.sigEvent)

        self.multiCh.SetNameTitle('MultiGraph_ch', 'MultiGraph_ch')
        self.multiCh.Add(self.adcCh)
        self.multiCh.Add(self.pedCh)
        self.multiCh.Add(self.cmnCh)
        self.multiCh.Add(self.sigCh)

    def BranchToVector(self, braName='DiaADC'):
        leng = self.selTree.Draw(braName,'','goff')
        if leng > 1000000:
            self.selTree.SetEstimate(leng)
            leng = self.selTree.Draw(braName,'','goff')
        vecSize = self.selTree.GetBranch(braName).GetLeaf(braName).GetLen()
        vector = [[self.selTree.GetV1()[ev*vecSize + ch] for ch in xrange(vecSize)] for ev in xrange(int(leng/vecSize))] if vecSize != 1 else [self.selTree.GetV1()[ev] for ev in xrange(int(leng))]
        return deepcopy(vector)


    def Get1DHistos(self, bat=False):
        if self.fidcut == 0 or self.fidcutCMN == 0:
            self.CreateFidCut()
        self.noiseHisto1D = self.noiseHisto2D.ProjectionY('{n}_all_Chs_Region_{r}'.format(n=self.noiseHisto2D.GetName(), r=self.nameChSelection), int(self.noiseHisto2D.GetXaxis().GetBinLowEdge(self.diamondNoiseChannels['ch_low'])), int(self.noiseHisto2D.GetXaxis().GetBinUpEdge(self.diamondNoiseChannels['ch_high'])+1))
        self.noiseHisto1D.GetXaxis().SetRangeUser(-32, 32)
        self.noiseHisto1D.SetFillColor(38)
        if bat: gROOT.SetBatch(True)
        canvas_name = 'c_noise_1D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('n')
        canvas.cd()
        self.noiseHisto1D.SetLineWidth(3)
        self.noiseHisto1D.GetYaxis().SetTitle('entries')
        self.FitHistogramGaus(self.noiseHisto1D, self.fitNoise, '{n}_fit'.format(n=self.noiseHisto1D.GetName()))
        gStyle.SetOptFit(0111)
        self.noiseHisto1D.Draw()
        name = 'histo1D_noise_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
        self.SaveCanvas(canvas, name)
        self.bla.append(canvas)
        if bat: gROOT.SetBatch(False)

        self.noiseCMNHisto1D = self.noiseCMNHisto2D.ProjectionY('{n}_all_Chs_Region_{r}'.format(n=self.noiseCMNHisto2D.GetName(), r=self.nameChSelection),
                                                          int(self.noiseCMNHisto2D.GetXaxis().GetBinLowEdge(
                                                              self.diamondNoiseChannels['ch_low'])),
                                                          int(self.noiseCMNHisto2D.GetXaxis().GetBinUpEdge(
                                                              self.diamondNoiseChannels['ch_high'])+1))
        self.noiseCMNHisto1D.GetXaxis().SetRangeUser(-32, 32)
        self.noiseCMNHisto1D.SetFillColor(38)
        if bat: gROOT.SetBatch(True)
        canvas_nameCMN = 'c_noiseCMNC_1D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                                                               l=self.diamondNoiseChannels['ch_low'],
                                                               h=self.diamondNoiseChannels['ch_high'])
        canvasCMN = self.CreateCanvas(canvas_nameCMN)
        gStyle.SetOptStat('n')
        canvasCMN.cd()
        self.noiseCMNHisto1D.SetLineWidth(3)
        self.noiseCMNHisto1D.GetYaxis().SetTitle('entries')
        self.FitHistogramGaus(self.noiseCMNHisto1D, self.fitNoiseCMNC, '{n}_fit'.format(n=self.noiseCMNHisto1D.GetName()))
        gStyle.SetOptFit(0111)
        self.noiseCMNHisto1D.Draw()
        nameCMN = 'histo1D_noiseCMNC_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
        self.SaveCanvas(canvasCMN, nameCMN)
        self.bla.append(canvasCMN)
        if bat: gROOT.SetBatch(False)

    def FitHistogramGaus(self, histo, fit, fit_name):
        mean = histo.GetMean()
        rms = histo.GetRMS()
        low = mean - 2*rms
        up = mean + 2*rms
        fit = TF1(fit_name, 'gaus', low, up)
        histo.Fit(fit_name, 'rq')
        fit.SetLineColor(kRed)
        fit.SetLineWidth(3)

    def Get2DMap(self, bat=False):
        if self.fidcut == 0 or self.fidcutCMN == 0:
            self.CreateFidCut()
        self.noiseHisto2D.GetXaxis().SetRangeUser(0, 128)
        self.noiseHisto2D.GetYaxis().SetRangeUser(-32, 32)
        self.noiseCMNHisto2D.GetXaxis().SetRangeUser(0, 128)
        self.noiseCMNHisto2D.GetYaxis().SetRangeUser(-32, 32)

        if bat: gROOT.SetBatch(True)
        canvas_name = 'c_noise_2D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('ne')
        canvas.cd()
        self.noiseHisto2D.Draw('colz')
        self.fidcut.Draw('same')
        self.bla.append(canvas)
        name = 'histo2D_noise_nonhit_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvas, name)

        canvasCMN_name = 'c_noiseCMNC_2D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        canvasCMN = self.CreateCanvas(canvasCMN_name)
        gStyle.SetOptStat('ne')
        canvasCMN.cd()
        self.noiseCMNHisto2D.Draw('colz')
        self.fidcutCMN.Draw('same')
        self.bla.append(canvasCMN)
        nameCMN = 'histo2D_noiseCMNC_nonhit_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvasCMN, nameCMN)
        if bat: gROOT.SetBatch(False)

    def Get2DMapFiducial(self, bat=False):
        if self.fidcut == 0 or self.fidcutCMN == 0:
            self.CreateFidCut()
        self.noiseHisto2D.GetXaxis().SetRangeUser(0, 128)
        self.noiseHisto2D.GetYaxis().SetRangeUser(-32, 32)
        self.noiseCMNHisto2D.GetXaxis().SetRangeUser(0, 128)
        self.noiseCMNHisto2D.GetYaxis().SetRangeUser(-32, 32)

        if bat: gROOT.SetBatch(True)
        canvas_name = 'c_noise_2D_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('ne')
        canvas.cd()
        self.noiseHisto2D.Draw('colz')
        self.fidcut.Draw('same')
        self.bla.append(canvas)
        name = 'histo2D_noise_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvas, name)

        canvasCMN_name = 'c_noiseCMNC_2D_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        canvasCMN = self.CreateCanvas(canvasCMN_name)
        gStyle.SetOptStat('ne')
        canvasCMN.cd()
        self.noiseCMNHisto2D.Draw('colz')
        self.fidcutCMN.Draw('same')
        self.bla.append(canvasCMN)
        nameCMN = 'histo2D_noiseCMNC_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvasCMN, nameCMN)
        if bat: gROOT.SetBatch(False)

        self.noiseHisto2D_fid = self.noiseHisto2D.Clone('{n}_region'.format(n=self.noiseHisto2D.GetName()))
        self.noiseHisto2D_fid.SetTitle(self.noiseHisto2D_fid.GetName())
        nbinsMap = (self.noiseHisto2D.GetNbinsY() * self.noiseHisto2D.GetNbinsX() + 1)
        for bin in xrange(1, nbinsMap):
            x, y, z = Long(0), Long(0), Long(0)
            if self.noiseHisto2D.GetBinContent(bin) > 0:
                self.noiseHisto2D.GetBinXYZ(bin, x, y, z)
                point = {'x': self.noiseHisto2D.GetXaxis().GetBinCenter(x), 'y': self.noiseHisto2D.GetYaxis().GetBinCenter(y)}
                # if not self.WindingIsPointInPoly(point):
                if not bool(int(self.fidcut.IsInside(point['x'], point['y']))):
                    self.noiseHisto2D_fid.SetBinContent(bin, 0)
        canvas_name = 'c_noise_2D_region_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                                                             l=int(self.diamondNoiseChannels['ch_low']),
                                                             h=int(self.diamondNoiseChannels['ch_high']))
        if bat: gROOT.SetBatch(True)
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('n')
        canvas.cd()
        self.noiseHisto2D_fid.Draw('colz')
        self.fidcut.Draw('same')
        self.bla.append(canvas)
        name = 'histo2D_noise_region_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                                                                l=int(self.diamondNoiseChannels['ch_low']),
                                                                h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvas, name)
        if bat: gROOT.SetBatch(False)

        self.noiseCMNHisto2D_fid = self.noiseCMNHisto2D.Clone('{n}_region'.format(n=self.noiseCMNHisto2D.GetName()))
        self.noiseCMNHisto2D_fid.SetTitle(self.noiseCMNHisto2D_fid.GetName())
        nbinsMap = (self.noiseCMNHisto2D.GetNbinsY() * self.noiseCMNHisto2D.GetNbinsX() + 1)
        for bin in xrange(1, nbinsMap):
            x, y, z = Long(0), Long(0), Long(0)
            if self.noiseCMNHisto2D.GetBinContent(bin) > 0:
                self.noiseCMNHisto2D.GetBinXYZ(bin, x, y, z)
                point = {'x': self.noiseCMNHisto2D.GetXaxis().GetBinCenter(x),
                         'y': self.noiseCMNHisto2D.GetYaxis().GetBinCenter(y)}
                # if not self.WindingIsPointInPoly(point):
                if not bool(int(self.fidcutCMN.IsInside(point['x'], point['y']))):
                    self.noiseCMNHisto2D_fid.SetBinContent(bin, 0)
        canvas_nameCMN = 'c_noiseCMNC_2D_region_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                                                                    l=int(self.diamondNoiseChannels['ch_low']),
                                                                    h=int(self.diamondNoiseChannels['ch_high']))
        if bat: gROOT.SetBatch(True)
        canvasCMN = self.CreateCanvas(canvas_nameCMN)
        gStyle.SetOptStat('n')
        canvasCMN.cd()
        self.noiseCMNHisto2D_fid.Draw('colz')
        self.fidcutCMN.Draw('same')
        self.bla.append(canvasCMN)
        nameCMN = 'histo2D_noiseCMNC_region_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                                                                       l=int(self.diamondNoiseChannels['ch_low']),
                                                                       h=int(self.diamondNoiseChannels['ch_high']))
        self.SaveCanvas(canvasCMN, nameCMN)
        if bat: gROOT.SetBatch(False)

    def CreateFidCut(self):
        if self.fidcut == 0 and self.fidcutCMN == 0:
            self.fidcut = TCutG('fidcut', 5)
            self.fidcut.SetVarX(self.noiseHisto2D.GetXaxis().GetName())
            self.fidcut.SetVarY(self.noiseHisto2D.GetYaxis().GetName())
            self.fidcut.SetPoint(0, self.diamondNoiseChannels['ch_low'], -32)
            self.fidcut.SetPoint(1, self.diamondNoiseChannels['ch_low'], 32)
            self.fidcut.SetPoint(2, self.diamondNoiseChannels['ch_high']+1, 32)
            self.fidcut.SetPoint(3, self.diamondNoiseChannels['ch_high']+1, -32)
            self.fidcut.SetPoint(4, self.diamondNoiseChannels['ch_low'], -32)
            self.fidcutCMN = TCutG('fidcutCMN', 5)
            self.fidcutCMN.SetVarX(self.noiseCMNHisto2D.GetXaxis().GetName())
            self.fidcutCMN.SetVarY(self.noiseCMNHisto2D.GetYaxis().GetName())
            self.fidcutCMN.SetPoint(0, self.diamondNoiseChannels['ch_low'], -32)
            self.fidcutCMN.SetPoint(1, self.diamondNoiseChannels['ch_low'], 32)
            self.fidcutCMN.SetPoint(2, self.diamondNoiseChannels['ch_high']+1, 32)
            self.fidcutCMN.SetPoint(3, self.diamondNoiseChannels['ch_high']+1, -32)
            self.fidcutCMN.SetPoint(4, self.diamondNoiseChannels['ch_low'], -32)
        if self.fidcut != 0 and self.fidcutCMN != 0:
            self.fidcut.SetLineColor(kRed)
            self.fidcut.SetLineWidth(3)
            self.fidcutCMN.SetLineColor(kRed)
            self.fidcutCMN.SetLineWidth(3)

    def SetDefaultChannels(self):
        self.diamondNoiseChannels = {'ch_low': 74, 'ch_high': 80}
        self.nameChSelection = 'Strip'

    def SetChannels(self):
        self.nameChSelection = raw_input('Enter the name for the channel selection: ')
        low, high = '', ''
        while True:
            if low == '':
                low = raw_input('Enter the lowest channel to analyse (int): ')
            if not self.CheckStringInt(low):
                print 'The channel entered is not an integer. Try again...'
                low = raw_input('Enter the lowest channel to analyse (int): ')
            elif 0 <= int(low) <= 127:
                low = int(low)
                break
            else:
                print 'The channel must be between [0, 127]. Try again...'
                low = raw_input('Enter the lowest channel to analyse (int): ')
        while True:
            if high == '':
                high = raw_input('Enter the highest channel to analyse (int): ')
            if not self.CheckStringInt(high):
                print 'The channel entered is not an integer. Try again...'
                high = raw_input('Enter the highest channel to analyse (int): ')
            elif low <= int(high) <= 127:
                high = int(high)
                break
            else:
                print 'The channel must be between [{l}, 127]. Try again...'.format(l=low)
                high = raw_input('Enter the highest channel to analyse (int): ')
        self.diamondNoiseChannels = {'ch_low': low, 'ch_high': high}
        self.CreateFidCut()

    def SetChannelsBatch(self, name, low, high):
        self.nameChSelection = str(name)
        self.diamondNoiseChannels = {'ch_low': low, 'ch_high': high}
        self.CreateFidCut()

    def CheckStringInt(self, s):
        if s[0] in ('-', '+'):
            return s[1:].isdigit()
        return s.isdigit()

    def CreateCanvas(self, name='c0', w=1024, h=1024):
        c = TCanvas(name, name, w, h)
        c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
        return deepcopy(c)

    def SaveCanvas(self, canvas, name):
        canvas.SaveAs(self.outputDir + '/root/{n}.root'.format(n=name))
        canvas.SaveAs(self.outputDir + '/Plots/{n}.png'.format(n=name))

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
    parser.add_option('-s', '--sourceDir', dest='source', default='./', type='string', help='source folder containing processed data of different runs')
    parser.add_option('-o', '--outputDir', dest='output', default='', type='string', help='output folder containing the analysed results')
    parser.add_option('-l', '--low', dest='low', default=1600, type='int', help='lower ADC to analyse (e.g. 1500)')
    parser.add_option('-t', '--top', dest='top', default=1665, type='int', help='topmost ADC to analyse (e.g. 1700)')

    (options, args) = parser.parse_args()
    run = int(options.run)
    source = str(options.source)
    output = str(options.output)
    low = int(options.low)
    top = int(options.top)

    z = SecondPulseID(run, source, output, low, top)
