import os, shutil, sys, progressbar
import numpy as np
import ROOT as ro
import ipdb

def IsInt(i):
	try:
		int(i)
		return True
	except ValueError:
		return False

def IsFloat(f):
	try:
		float(f)
		return True
	except ValueError:
		return False

def RoundInt(n, nptype='int32'):
	val = np.floor(np.add(n, 0.5, dtype='f8'), dtype='f8').astype(nptype)
	if nptype.lower().startswith('i'):
		return int(val)
	elif nptype.lower().startswith('f'):
		return float(val)
	return val

def IsOnlyInt(f):
	if IsFloat(f):
		return RoundInt(f) == f
	else:
		return False

def CreateDefaultSettingsFile(diri, run_no, events, dut_name='default', dut_volt=0, ev_ini=0, num_evs_ana=0, dia_input=0, dia_sat=4095, max_trans_clust=10, num_highest_trans=5, chi2=5):
	if not os.path.isdir(diri):
		os.makedirs(diri)
	with open(diri + '/settings.{r}.ini'.format(r=run_no), 'w') as f:
		f.write('runNo = {r}\n'.format(r=run_no))
		f.write('Events = {e}\n'.format(e=events))
		f.write('event_start = {ei}\n'.format(ei=ev_ini))
		f.write('Events_ana = {ea}\n'.format(ea=events-ev_ini if num_evs_ana == 0 else num_evs_ana))
		f.write('repeaterCardNo = 2\n')
		f.write('voltage = {v}\n'.format(v=dut_volt))
		f.write('diamondName = {dn}\n'.format(dn=dut_name))
		f.write('currentBegin = 0\n')
		f.write('currentEnd = 0\n\n')
		f.write('diamondPattern = {0,50,0,127}\n\n')
		f.write('Iter_Size = 500\n\n')
		f.write('dia_input = {d}\n\n'.format(d=dia_input))
		f.write('dia_saturation = {d}\n\n'.format(d=dia_sat))
		f.write('Dia_channel_screen_channels = {}\n\n')
		f.write('Dia_channel_not_connected = {}\n\n')
		f.write('Dia_channel_noisy = {}\n\n')
		f.write('D0X_channel_screen_channels = {0,255}\n')
		f.write('D0Y_channel_screen_channels = {0,1,2,127,255}\n')
		f.write('D1X_channel_screen_channels = {0,1,214,215,232,233}\n')
		f.write('D1Y_channel_screen_channels = {127,255}\n')
		f.write('D2X_channel_screen_channels = {0,255}\n')
		f.write('D2Y_channel_screen_channels = {0,127,255}\n')
		f.write('D3X_channel_screen_channels = {0,1,2,3,127,255}\n')
		f.write('D3Y_channel_screen_channels = {0,1,2,3,246,247,253,254,255}\n\n')
		f.write('Si_Pedestal_Hit_Factor = 5\n')
		f.write('Di_Pedestal_Hit_Factor = 3\n\n')
		f.write('DO_CMC = 1\n\n')
		f.write('clusterSeedFactors = {16,20,20,20,18,18,14,10,5}\n')
		f.write('clusterHitFactors = {10,16,14,14,14,14,12,10,3}\n\n')
		f.write('UseAutoFidCut = 0\n')
		f.write('nDiamonds = 1\n')
		f.write('eta_corr_limit = 0.0005\n')
		f.write('si_avg_fidcut_xlow = 0\n')
		f.write('si_avg_fidcut_xhigh = 255\n')
		f.write('si_avg_fidcut_ylow = 0\n')
		f.write('si_avg_fidcut_yhigh = 255\n')
		f.write('selectionFidCut = {0-255,0-255}\n')
		f.write('RerunSelection = 0\n\n')
		f.write('DetectorsToAlign = 2\n')
		f.write('alignment_training_fidcuts = {1}\n')
		f.write('alignment_training_track_number = {e}\n'.format(e=int(round(num_evs_ana)/10.0)))
		f.write('alignment_training_method = 1\n')
		f.write('AlignmentIgnoreChannelDia = {0,127}\n')
		f.write('alignment_chi2 = {ch}\n\n'.format(ch=chi2))
		f.write('Double_t detectorD0Z = 0.725\n')
		f.write('Double_t detectorD1Z = 1.625\n')
		f.write('Double_t detectorD2Z = 18.725\n')
		f.write('Double_t detectorD3Z = 19.625\n')
		f.write('Double_t detectorDiaZ = 10.2\n\n')
		f.write('TransparentAlignment = 0\n\n')
		f.write('transparentChi2 = 0\n')
		f.write('yOffset3D = 0\n')
		f.write('max_transparent_cluster_size = {d}\n'.format(d=max_trans_clust))
		f.write('num_highest_transparent_cluster = {d}\n\n'.format(d=num_highest_trans))
		f.write('pulse_height_di_max = 4000\n')
		f.write('pulse_height_num_bins = 200\n')
		f.write('3dShortAnalysis = 0;\n')
		f.write('3dLongAnalysis = 0;\n')
		f.write('3dTransparentAnalysis = 0;\n')

def CreateDirectoryIfNecessary(dir):
	if not os.path.isdir(dir):
		print 'Creating directory', dir, '...', ; sys.stdout.flush()
		os.makedirs(dir)
		print 'Done'

def Set_voltage(s_file):
	Replace_Settings_Line(s_file, 'voltage')

def Mark_channels_as_NC(s_file):
	Replace_Settings_Line(s_file, 'Dia_channel_not_connected')

def Mark_channels_as_Screened(s_file):
	Replace_Settings_Line(s_file, 'Dia_channel_screen_channels')

def Mark_channels_as_Noisy(s_file):
	Replace_Settings_Line(s_file, 'Dia_channel_noisy')

def Set_Diamond_pattern(s_file):
	Replace_Settings_Line(s_file, 'diamondPattern')

def Set_Diamond_name(s_file):
	Replace_Settings_Line(s_file, 'diamondName')

def Mark_channels_as_no_alignment(s_file):
	with open(s_file, 'r') as f0:
		for line in f0:
			if line.lower().startswith('Dia_channel_not_connected'.lower()):
				print 'The current value(s) for not connected channels are:'
				print line.split('=')[1].split('\n')[0]
			if line.lower().startswith('Dia_channel_screen_channels'.lower()):
				print 'The current value(s) for screened channels are:'
				print line.split('=')[1].split('\n')[0]
			if line.lower().startswith('Dia_channel_noisy'.lower()):
				print 'The current value(s) for noisy channels are:'
				print line.split('=')[1].split('\n')[0]
	Replace_Settings_Line(s_file, 'AlignmentIgnoreChannelDia')

def Select_fiducial_region(s_file):
	xlow = Replace_Settings_Line(s_file, 'si_avg_fidcut_xlow')
	xhigh = Replace_Settings_Line(s_file, 'si_avg_fidcut_xhigh')
	ylow = Replace_Settings_Line(s_file, 'si_avg_fidcut_ylow')
	yhigh = Replace_Settings_Line(s_file, 'si_avg_fidcut_yhigh')
	num_patt = 0
	with open(s_file, 'r') as f0:
		for line in f0:
			if line.lower().startswith('diamondPattern'.lower()):
				num_patt += 1
	if num_patt == 1:
		Replace_Settings_Line(s_file, 'selectionFidCut', action='value', value=('{' + xlow + '-' + xhigh + ',' + ylow + '-' + yhigh + '}').replace(' ', ''))

def Replace_Settings_Line(s_file, setting_option, action='user', value=''):
	if not os.path.isfile(s_file):
		print 'The file', s_file, 'does not exist'
		return
	updated = False
	with open(s_file, 'r') as f0:
		with open('temp00.ini', 'w') as f1:
			for line in f0:
				if not line.lower().startswith(setting_option.lower()):
					f1.write(line)
				else:
					if action == 'user':
						new_values = Get_From_User_Line(setting_option, line.split('=')[1].split('\n')[0])
					elif action == 'even':
						new_values = Only_Even_Channels(line)
					elif action == 'odd':
						new_values = Only_Odd_Channels(line)
					elif action == 'value':
						new_values = value
					f1.write(setting_option + ' = ' + new_values + '\n')
					updated = True
	if not updated:
		with open('temp00.ini', 'a') as f1:
			if action == 'user':
				new_values = Get_From_User_Line(setting_option, update=False)
			elif action == 'even':
				new_values = Only_Even_Channels('Dia_channel_screen_channels = {1,127}')
			elif action == 'odd':
				new_values = Only_Odd_Channels('Dia_channel_screen_channels = {0,126}')
			f1.write(setting_option + ' = ' + new_values + '\n')
	os.remove(s_file)
	shutil.move('temp00.ini', s_file)
	if setting_option.startswith('si_avg_fidcut'):
		return new_values

def Get_From_User_Line(option, old_value='', update=True):
	if update:
		print 'The current value(s) for:', option, 'is(are):'
		print old_value
		temp = raw_input('Enter the new value(s) for ' + option + ' in the same format as above (Press enter to leave it as it is): ')
		if temp == '' or temp == '\n':
			temp = old_value
	else:
		temp = raw_input('Enter the value(s) for ' + option + ' in the correct format: ')
	return temp.replace(' ', '')

def Get_From_User_Value(option, type='float', old_value=0, update=True):
	if update:
		print 'The current value for:', option, 'is:'
		print '', old_value
		temp = raw_input('Enter the new value for ' + option + ' which should be a ' + type + ' (Press enter to leave it as it is): ')
		temp = old_value if temp in ['', '\n'] else int(temp) if type.lower().startswith('int') and IsInt(temp) else float(temp) if type.lower().startswith('float') and IsFloat(temp) else 0
	else:
		temp = raw_input('Enter the value for ' + option + ' which should be a ' + type + ': ')
		temp = int(temp) if type.lower().startswith('int') and IsInt(temp) else float(temp) if type.lower().startswith('float') and IsFloat(temp) else 0
	return temp

def Only_Even_Channels(old_value):
	channel_str = old_value[old_value.find('{') + 1:old_value.find('}')].split(',')
	channel_str = ['0', '126'] if len(channel_str) == 1 and channel_str[0] == '' else channel_str
	return Modify_String_Even_or_Odd(channel_str, 'even')
	
def Only_Odd_Channels(old_value):
	channel_str = old_value[old_value.find('{') + 1:old_value.find('}')].split(',')
	channel_str = ['1', '127'] if len(channel_str) == 1 and channel_str[0] == '' else channel_str
	return Modify_String_Even_or_Odd(channel_str, 'odd')
	
def Modify_String_Even_or_Odd(channel_old, type='even'):
	if type == 'even':
		modulo = 1
		if '1' != channel_old[0] and not '1-' in channel_old[0]:
			if '0' != channel_old[0] and not '0-' in channel_old[0]:
				channel_old.insert(0, '1')
		if '127' != channel_old[-1] and not '-127' in channel_old[-1]:
			channel_old.append('127')
		channel_old.append('129')
	else:
		modulo = 0
		if '0' != channel_old[0] and not '0-' in channel_old[0]:
			channel_old.insert(0, '0')
		if '126' != channel_old[-1] and not '-126' in channel_old[-1]:
			if '127' != channel_old[-1] and not '-127' in channel_old[-1]:
				channel_old.append('126')
		channel_old.append('128')
	channel_new = ''
	for i in xrange(1, len(channel_old)):
		if IsInt(channel_old[i - 1]):
			prev = int(channel_old[i - 1])
			if IsInt(channel_old[i]):
				th = int(channel_old[i])
				if th - prev < 2:
					channel_new += str(prev) + ','
				elif prev % 2 == modulo:
					for ch in xrange(prev, th, 2):
						channel_new += str(ch) + ','
				else:
					channel_new += str(prev) + ','
					for ch in xrange(prev + 1, th, 2):
						channel_new += str(ch) + ','
			else:
				temp = channel_old[i].split('-')
				if IsInt(temp[0]):
					th = int(temp[0])
					if th - prev < 2:
						channel_new += str(prev) + ','
					elif prev % 2 == modulo:
						for ch in xrange(prev, th, 2):
							channel_new += str(ch) + ','
					else:
						channel_new += str(prev) + ','
						for ch in xrange(prev + 1, th, 2):
							channel_new += str(ch) + ','
		else:
			channel_new += channel_old[i - 1] + ','
			prev = int(channel_old[i - 1].split('-')[-1])
			if not IsInt(channel_old[i]):
				th = int(channel_old[i].split('-')[0])
			else:
				th = int(channel_old[i])
			if th - prev >= 2:
				if prev % 2 == modulo:
					for ch in xrange(prev + 2, th, 2):
						channel_new += str(ch) + ','
				else:
					for ch in xrange(prev + 1, th, 2):
						channel_new += str(ch) + ','

	channel_new = channel_new[:-1]
	return '{' + channel_new + '}'

def RecreateLink(source, dest, name, doSymlink=True, doCopy=True):
	successful = True
	if os.path.isdir(source):
		if doSymlink:
			successful = RecreateSoftLink(source, dest, name, type='dir', doCopy=doCopy)
		elif doCopy:
			if os.path.islink(dest + '/' + name):
				os.unlink(dest + '/' + name)
			elif os.path.isdir(dest + '/' + name):
				shutil.rmtree(dest + '/' + name, True)
			shutil.copytree(source, dest + '/' + name)
			successful = True
		else:
			successful = False
			print 'Could not link', source, 'to', dest, 'and copy was disabled'

	elif os.path.isfile(source):
		if os.path.isfile(dest + '/' + name):
			os.unlink(dest + '/' + name)
		if doSymlink:
			successful = RecreateSoftLink(source, dest, name, type='file', doCopy=doCopy)
		if not doSymlink or not successful:
			if os.path.isfile(dest + '/' + name):
				os.unlink(dest + '/' + name)
			try:
				os.link(source, dest + '/' + name)
				successful = True
			except OSError:
				if doCopy:
					shutil.copy2(source, dest + '/' + name)
					successful = True
				else:
					successful = False
					print 'Could not link', source, 'to', dest, 'and copy was disabled'
	else:
		successful = False
	return successful

def RecreateSoftLink(source, dest, name, type='dir', doCopy=True):
	if type == 'dir':
		if os.path.isdir(source):
			if os.path.islink(dest + '/' + name):
				os.unlink(dest + '/' + name)
			elif os.path.isdir(dest + '/' + name):
				shutil.rmtree(dest + '/' + name, True)
			try:
				os.symlink(source, dest + '/' + name)
				return True
			except OSError:
				if doCopy:
					shutil.copytree(source, dest + '/' + name)
					return True
	else:
		if os.path.isfile(source):
			if os.path.isfile(dest + '/' + name):
				os.unlink(dest + '/' + name)
			try:
				os.symlink(source, dest + '/' + name)
				return True
			except OSError:
				if doCopy:
					shutil.copy2(source, dest + '/' + name)
					return True
	return False

def CloseSubprocess(p, stdin=False, stdout=False):
	pid = p.pid
	if stdin:
		p.stdin.close()
	if stdout:
		p.stdout.close()
	if p.wait() is None:
		print 'Could not terminate subprocess... forcing termination'
		p.kill()
		if p.wait() is None:
			ExitMessage('Could not kill subprocess... quitting')
	try:
		os.kill(pid, 0)
	except OSError:
		pass
	else:
		print 'The subprocess is still running. Killing it with os.kill'
		os.kill(pid, 15)
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			ExitMessage('The process does not die... quitting program')
	del p, pid
	p = None
	print 'Subprocess terminated successfully'

def ReturnRGB(val, min, max):
		hue = (val - min) * 1.0 / float(max - min) if max != min else 0
		huep = hue * 6.0
		i = int(huep)
		v1 = 0.
		v2 = 1. - huep + i
		v3 = huep - i
		if i == 0:
			r, g, b = 1., v3, v1
		elif i == 1:
			r, g, b = v2, 1., v1
		elif i == 2:
			r, g, b = v1, 1., v3
		elif i == 3:
			r, g, b = v1, v2, 1.
		elif i == 4:
			r, g, b = v3, v1, 1.
		else:
			r, g, b = 1., v1, v2
		return r, g, b

def Correct_Path(path, times=2):
	abs_path = ''
	if path[0] == '~':
		abs_path += os.path.expanduser('~')
		abs_path += path[1:]
	elif os.path.isabs(path):
		abs_path += path
	else:
		abs_path += os.path.abspath(path)
	if times != 1:
		return Correct_Path(abs_path, 1)
	return abs_path

def DeleteDirectoryContents(dir):
	if os.path.isdir(dir):
		print 'Deleting old analysis...', ; sys.stdout.flush()
		for file in os.listdir(dir):
			file_path = os.path.join(dir, file)
			try:
				if os.path.islink(file_path) or os.path.isfile(file_path):
					os.unlink(file_path)
				elif os.path.isdir(file_path):
					shutil.rmtree(file_path)
			except Exception as e:
				print(e)
		print 'Done'

def AddLineToStats(canvas, key, value=0, samplelinekey='Mean'):
	if canvas:
		ps = canvas.GetPrimitive('stats')
		ps.SetName('mystats')
		lol = ps.GetListOfLines()
		sampleline = ps.GetLineWith(samplelinekey)
		line = ro.TLatex(0, 0, '{k} = {v:.2f}'.format(k=key, v=value))
		line.SetTextAlign(sampleline.GetTextAlign())
		line.SetTextAngle(sampleline.GetTextAngle())
		line.SetTextColor(sampleline.GetTextColor())
		line.SetTextFont(sampleline.GetTextFont())
		line.SetTextSize(sampleline.GetTextSize())
		lol.Add(line)
		canvas.Modified()

def SetX1X2NDC(histo, x1, x2, statsname='stats'):
	if histo.FindObject(statsname):
		histo.FindObject(statsname).SetX1NDC(x1)
		histo.FindObject(statsname).SetX2NDC(x2)
		ro.gPad.Update()

def SetLegendX1X2Y1Y2(legend, x1, x2, y1, y2):
	legend.SetX1NDC(x1)
	legend.SetX2NDC(x2)
	legend.SetY1NDC(y1)
	legend.SetY2NDC(y2)
	ro.gPad.Update()

def CreateProgressBarUtils(maxVal=1):
	widgets = [
		'Processed: ', progressbar.Counter(),
		' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
		' ', progressbar.Bar(marker='>'),
		' ', progressbar.Timer(),
		' ', progressbar.ETA()
		# ' ', progressbar.AdaptativeETA(),
		#  ' ', progressbar.AdaptativeTransferSpeed()
	]
	bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)
	return bar

def SetDefault2DStats(histo):
	if histo.FindObject('stats'):
		histo.FindObject('stats').SetOptStat(11)
		histo.FindObject('stats').SetX1NDC(0.7)
		histo.FindObject('stats').SetX2NDC(0.9)
		histo.FindObject('stats').SetY1NDC(0.9)
		histo.FindObject('stats').SetY2NDC(0.975)
		ro.gPad.Update()
	if histo.FindObject('palette'):
		histo.FindObject('palette').SetX1NDC(0.87)
		histo.FindObject('palette').SetX2NDC(0.92)
		ro.gPad.Update()
		
def SetDefault1DStats(histo):
	if histo.FindObject('stats'):
		histo.FindObject('stats').SetOptStat(112211)
		histo.FindObject('stats').SetX1NDC(0.6)
		histo.FindObject('stats').SetX2NDC(0.9)
		histo.FindObject('stats').SetY1NDC(0.6)
		histo.FindObject('stats').SetY2NDC(0.9)
		ro.gPad.Update()

def FindLeafInTree(tree, leaf):
	if tree.FindLeaf(leaf):
		return True
	return False

def ExitMessage(txt, code=os.EX_SOFTWARE):
    print '##########'
    print txt
    print '##########'
    sys.exit(code)

