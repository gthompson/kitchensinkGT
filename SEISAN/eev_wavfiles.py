import sys
import StreamGT

def main():

	filenrlis = sys.argv[1]
	print filenrlis
	file=open(filenrlis,'r') 

	lines = file.readlines()

	for line in lines:
		wavfile = line[6:].strip()

		# read the wav file as a stream object
                st = StreamGT()
		st.read(wavfile)
		print wavfile
                st.statistics()               

		print raw_input('<ENTER> to move to next file')

if __name__ == "__main__":
    main()


