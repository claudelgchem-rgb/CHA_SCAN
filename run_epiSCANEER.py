import glob
import os
import shutil
import sys

from epiSCANEER.epiSCANEER import *
from epiSCANEER import msa, coe, iscalc

input_path = "./input/"
output_path = "./output/"

# Optional command line overrides: python run_epiSCANEER.py [input_dir] [output_dir]
if len(sys.argv) > 1:
	input_path = sys.argv[1]
if len(sys.argv) > 2:
	output_path = sys.argv[2]

if not os.path.isdir(output_path):
	os.makedirs(output_path)

file_list = sorted(glob.glob(os.path.join(input_path, "*.aln")))
if not file_list:
	sys.exit("No *.aln file found in %s" % os.path.abspath(input_path))

for msa_path in file_list:
	msa_file = os.path.basename(msa_path)
	prefix = msa_file.split('.')[0]

	base_pth = os.path.join(output_path, prefix)
	if not os.path.isdir(base_pth):
		os.makedirs(base_pth)
	elif os.path.isfile(os.path.join(base_pth, "%s.txt" % prefix)):
		print("already calculated! skip! (%s)" % prefix)
		continue

	# Work on a copy of the MSA so that every result file is written to output/
	work_msa_path = os.path.join(base_pth, "%s.aln" % prefix)
	if os.path.abspath(work_msa_path) != os.path.abspath(msa_path):
		shutil.copyfile(msa_path, work_msa_path)

	# Loading MSA
	pm = msa.ProcMsa("tmp", work_msa_path, "tmp", "tmp")
	pm.parse()
	msa_dic = build_msa(work_msa_path, prefix)

	# Construct co-evolutionary network and Calculating SCI
	calc_epiSCI(base_pth, msa_dic, pm, prefix)

	print("done: %s -> %s" % (prefix, os.path.join(base_pth, "%s.txt" % prefix)))
