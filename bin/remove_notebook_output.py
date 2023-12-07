import sys
from nbformat import read, write

def Remove_Output(Book):
    for cell in Book.cells:
        if hasattr(cell, "outputs"):
            cell.outputs = []
        if hasattr(cell, "prompt_number"):
            del cell["prompt_number"]

print(len(sys.argv))
for i, arg in enumerate(sys.argv):
    print(i, arg)

if len(sys.argv)<2:
    print(f"Usage: {sys.argv[0]} notebook_in [notebook_out]")
    exit()
elif len(sys.argv)==2:
    infile = sys.argv[1]
    outfile = infile
elif len(sys.argv)==3:
    infile = sys.argv[1]
    outfile = sys.argv[2]
Book= read(open(infile), 4)
Remove_Output(Book)
write(Book, open(outfile, "w"), 4)