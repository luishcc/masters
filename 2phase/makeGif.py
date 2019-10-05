import os, sys
import datetime
import imageio
from pprint import pprint
import time
import datetime

e = sys.exit
 
cwd = os.getcwd()

newdir = cwd + '/2D/result'
#newdir = cwd + '/3D/result'

os.chdir(newdir)

def create_gif(filenames, duration):
	images = []

	for filename in filenames:
		images.append(imageio.imread(filename))
	output_file = 'Gif-%s.gif' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')
	os.chdir(cwd)
	imageio.mimsave(output_file, images, duration=duration)

if __name__ == "__main__":
	script = sys.argv.pop(0)
	duration = 0.01
	filenames = sorted(filter(os.path.isfile, [x for x in os.listdir(newdir) if x.endswith(".jpg")]),
					   key=lambda p: os.path.exists(p) and os.stat(p).st_mtime or
									 time.mktime(datetime.now().timetuple()))

	create_gif(filenames, duration)
