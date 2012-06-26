
import cProfile,subprocess
#from scrolltest import Compressor
from recipsample import Compressor
command ="""Compressor()"""
cProfile.runctx(command,globals(),locals(),filename="profile.txt")
subprocess.check_call(['runsnake','profile.txt'])