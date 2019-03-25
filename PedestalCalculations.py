# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
# sys.path.append('/home/sandiego/rd42AnalysisTools')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
# from Settings import Settings
from Utils import *
from PedestalDeviceCalculations import PedestalDeviceCalculations
import time
import numpy as np
import subprocess as subp
# from threading import Thread
# from multiprocessing import Manager, sharedctypes
import multiprocessing as mp
import multiprocessing.sharedctypes as mpsc
import ctypes
import ipdb

__author__ = 'DA'

class PedestalCalculations:
    def __init__(self, origTree, outdir, run_no, sliding_length=50, hit_fact=3, seed_fact=5):
        print 'Creating PedestalCalculations instance'
        self.rootTree = origTree
        self.out_dir = outdir
        self.run_no = run_no
        self.events = self.rootTree.GetEntries()
        self.raw_dut_branch = 'DiaADC'
        self.chs = self.rootTree.GetLeaf(self.raw_dut_branch).GetLen()
        self.dut_np_data_type, self.dut_root_data_type = 'uint16', 's'
        self.hit_fact = hit_fact
        self.seed_fact = seed_fact



        # self.out_dir = self.settings.output_dir
        # self.sub_dir = self.settings.sub_dir
        # self.file_name = self.settings.file_name
        # self.run = self.settings.run
        # self.tree_name = self.settings.tree_name
        self.num_parallel = RoundInt(mp.cpu_count() / 2.0)
        self.slide_leng = sliding_length
        self.cm_cut = 4
        # self.ana_dir = self.settings.analysis_path
        self.ped_branches = ['diaHitChs', 'diaSeedChs', 'diaPedestalMean', 'diaPedestalSigma', 'diaCmChs', 'diaCm', 'diaHitChsCmc', 'diaSeedChsCmc', 'diaPedestalMeanCmc', 'diaPedestalSigmaCmc']
        self.devices = ['dut']
        self.device_ped_process = {}
        self.outFile = ro.TFile()
        self.outTree = ro.TTree()

        self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs = np.zeros(self.chs, '?'), np.zeros(self.chs, '?'), np.zeros(self.chs, '?')
        self.GetScreenedChannels()

        # self.dic_device_to_pos = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}

        # Temporary numpy arrays to create the ctype vectors for multiprocessing with shared memory

        dut_ADC_mean_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_sigma_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_is_ped_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_hit_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_seed_all = np.zeros((self.chs, self.events), dtype='uint8')

        # Temporary numpy arrays to create the ctype vectors for multiprocessing with shared memory for cmc calculations
        dut_ADC_chs_cm_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_cm_all = np.zeros(self.events, dtype='float32')
        dut_ADC_mean_cmc_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_sigma_cmc_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_is_ped_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_hit_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_seed_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')

        # The ctype vectors that will be shared in the multiprocessing
        self.dut_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_mean_all))
        self.dut_ADC_sigma_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_sigma_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_sigma_all))
        self.dut_ADC_is_ped_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_ped_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_ped_all))
        self.dut_ADC_is_hit_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_hit_all))
        self.dut_ADC_is_seed_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_seed_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_seed_all))
        del dut_ADC_mean_all, dut_ADC_sigma_all, dut_ADC_is_ped_all, dut_ADC_is_hit_all, dut_ADC_is_seed_all

        # The ctype vectors that will be shared in the multiprocessing for cmc calculations
        self.dut_ADC_chs_cm_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_chs_cm_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_chs_cm_all))
        self.dut_ADC_cm_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_cm_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_cm_all))
        self.dut_ADC_mean_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_mean_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_mean_cmc_all))
        self.dut_ADC_sigma_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_sigma_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_sigma_cmc_all))
        self.dut_ADC_is_ped_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_ped_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_ped_cmc_all))
        self.dut_ADC_is_hit_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_hit_cmc_all))
        self.dut_ADC_is_seed_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_seed_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_seed_cmc_all))
        del dut_ADC_chs_cm_all, dut_ADC_cm_all, dut_ADC_mean_cmc_all, dut_ADC_sigma_cmc_all, dut_ADC_is_ped_cmc_all, dut_ADC_is_hit_cmc_all, dut_ADC_is_seed_cmc_all

        # ctype vectors with all the ADC events for telescope and dut. The each process will read from this vectors
        self.dut_ADC_all_mp = self.LoadDutADCs()

        # numpy arrays to fill tree
        self.dut_ADC_is_ped = np.zeros(self.chs, '?')
        self.dut_ADC_is_hit = np.zeros(self.chs, '?')
        self.dut_ADC_is_seed = np.zeros(self.chs, '?')
        self.dut_ADC_mean = np.zeros(self.chs, 'float32')
        self.dut_ADC_sigma = np.zeros(self.chs, 'float32')
        self.dut_ADC_cm_chs = np.zeros(self.chs, '?')
        self.dut_ADC_cm = np.zeros(1, 'float32')
        self.dut_ADC_is_ped_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_is_hit_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_is_seed_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_mean_cmc = np.zeros(self.chs, 'float32')
        self.dut_ADC_sigma_cmc = np.zeros(self.chs, 'float32')

        # branches to fill tree
        self.b1, self.b2, self.b3, self.b4, self.b5, self.b6, self.b7, self.b8, self.b9, self.b10, self.b11, self.b12, self.b13, self.b14, self.b15, self.b16, self.b17 = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

    def GetScreenedChannels(self):
        list_bra = ['diaNcChs', 'diaNoisyChs', 'diaMaskedChs']
        Draw_Branches_For_Get_Val(self.rootTree, list_bra)
        Get_Branches_Value_To_Numpy(self.rootTree, list_bra, [self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs], 1, self.chs)

    def LoadDutADCs(self):
        # self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        # self.rootTree.SetEstimate(256 * self.rootTree.GetEntries())
        time0 = time.time()
        print 'Getting all events for DUT...', ; sys.stdout.flush()

        Draw_Branches_For_Get_Val(self.rootTree, [self.raw_dut_branch], 0, self.events, option='goff')
        temp = np.zeros((self.chs, self.events), dtype=self.dut_np_data_type)
        Get_Branches_Value_To_Numpy(self.rootTree, [self.raw_dut_branch], [temp], self.events, self.chs)
        print 'Done in', time.time() - time0, 'seconds'
        return mp.Array(np.ctypeslib.as_ctypes(temp)._type_, np.ctypeslib.as_ctypes(temp))

    def CalculateDevicesPedestals(self):
        jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
        for row in jobs_chunks:
            self.device_ped_process = []
            for it, device in enumerate(row):
                if it == len(row) - 1:
                    print 'Calculating pedestals for', device, '. The progress is shown below:'
                    temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)

                else:
                    print 'Calculating pedestals for', device
                    temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                temp.start()
                self.device_ped_process.append(temp)
            for j in self.device_ped_process:
                j.join()

        # self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='UPDATE')
        self.outFile = ro.TFile('{d}/{r}/pedestalData.{r}.{s}.root'.format(d=self.out_dir, r=self.run_no, s=self.slide_leng), 'RECREATE')
        self.outTree = ro.TTree('pedestalTree', 'pedestalTree')
        self.SetBranches()
        self.FillTree()


    def SetBranches(self):
        print 'Setting branches in tree...', ; sys.stdout.flush()
        self.b6 = self.outTree.Branch('diaPedChs', self.dut_ADC_is_ped, 'diaPedChs[{chs}]/O'.format(chs=self.chs))
        self.b7 = self.outTree.Branch('diaHitChs', self.dut_ADC_is_hit, 'diaHitChs[{chs}]/O'.format(chs=self.chs))
        self.b8 = self.outTree.Branch('diaSeedChs', self.dut_ADC_is_seed, 'diaSeedChs[{chs}]/O'.format(chs=self.chs))
        self.b9 = self.outTree.Branch('diaPedestalMean', self.dut_ADC_mean, 'diaPedestalMean[{f}]/F'.format(f=self.chs))
        self.b10 = self.outTree.Branch('diaPedestalSigma', self.dut_ADC_sigma, 'diaPedestalSigma[{f}]/F'.format(f=self.chs))
        self.b11 = self.outTree.Branch('diaCmChs', self.dut_ADC_cm_chs, 'diaCmChs[{f}]/O'.format(f=self.chs))
        self.b12 = self.outTree.Branch('diaCm', self.dut_ADC_cm, 'diaCm/F')
        self.b13 = self.outTree.Branch('diaPedChsCmc', self.dut_ADC_is_ped_cmc, 'diaPedChsCmc[{chs}]/O'.format(chs=self.chs))
        self.b14 = self.outTree.Branch('diaHitChsCmc', self.dut_ADC_is_hit_cmc, 'diaHitChsCmc[{chs}]/O'.format(chs=self.chs))
        self.b15 = self.outTree.Branch('diaSeedChsCmc', self.dut_ADC_is_seed_cmc, 'diaSeedChsCmc[{chs}]/O'.format(chs=self.chs))
        self.b16 = self.outTree.Branch('diaPedestalMeanCmc', self.dut_ADC_mean_cmc, 'diaPedestalMeanCmc[{f}]/F'.format(f=self.chs))
        self.b17 = self.outTree.Branch('diaPedestalSigmaCmc', self.dut_ADC_sigma_cmc, 'diaPedestalSigmaCmc[{f}]/F'.format(f=self.chs))
        print 'Done'

    def FillTree(self):
        print 'Filling tree...'
        bar = CreateProgressBarUtils(self.events)
        bar.start()
        for ev in xrange(self.events):
            self.outTree.GetEntry(ev)
            self.LoadArrays(ev)
            bar.update(ev + 1)
        self.outTree.Write()
        bar.finish()
        self.outFile.Close()

    def LoadArrays(self, ev):
        np.putmask(self.dut_ADC_is_ped, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_ped_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_is_hit, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_hit_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_is_seed, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_seed_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_mean, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_mean_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_sigma, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_sigma_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_cm_chs, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_chs_cm_all_mp.get_obj())[:, ev].astype('?'))
        self.dut_ADC_cm.fill(np.ctypeslib.as_array(self.dut_ADC_cm_all_mp.get_obj())[ev])
        np.putmask(self.dut_ADC_is_ped_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_ped_cmc_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_is_hit_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_hit_cmc_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_is_seed_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_seed_cmc_all_mp.get_obj())[:, ev].astype('?'))
        np.putmask(self.dut_ADC_mean_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_mean_cmc_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_sigma_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_sigma_cmc_all_mp.get_obj())[:, ev])
        for branch in [self.b6, self.b7, self.b8, self.b9, self.b10, self.b11, self.b12, self.b13, self.b14, self.b15, self.b16, self.b17]:
            branch.Fill()


def main():
    z = PedestalCalculations()


if __name__ == '__main__':
    main()
