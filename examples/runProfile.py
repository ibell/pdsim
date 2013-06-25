
import cProfile,subprocess
from scrolltest import Compressor
#from recipsample import Compressor
#from PUrecip import Compressor
Compressor(OneCycle = False)
#cProfile.runctx("""Compressor(OneCycle=False)""",globals(),locals(),filename="profile.txt")
#subprocess.check_call(['runsnake','profile.txt'])