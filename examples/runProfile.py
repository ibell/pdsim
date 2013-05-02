
import cProfile,subprocess
from scrolltest import Compressor
#from recipsample import Compressor
#from PUrecip import Compressor
cProfile.runctx("""Compressor(OneCycle=True)""",globals(),locals(),filename="profile.txt")
subprocess.check_call(['runsnake','profile.txt'])