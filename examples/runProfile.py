
import cProfile,subprocess
from scrolltest import Compressor
#from recipsample import Compressor
#from PUrecip import Compressor
cProfile.runctx("""Compressor(TTSE=False)""",globals(),locals(),filename="profile.txt")
cProfile.runctx("""Compressor(TTSE=True)""",globals(),locals(),filename="profile2.txt")
subprocess.check_call(['runsnake','profile.txt'])
subprocess.check_call(['runsnake','profile2.txt'])