import numpy
import csv
import argparse
import sys
from PyQt4 import QtGui
from argparseui import ArgparseUi

def main():
	parser = argparse.ArgumentParser(description="Extract information from NeuronScience files")
	parser.add_argument("-swc", required=True, help="Path to the SWC file containing dendrite data.")
	parser.add_argument("-spine", required=True, help="Path to the spine TXT file.")
	parser.add_argument("-output", default="output.csv", help="Path to the output CSV that will be written.")

	app = QtGui.QApplication(sys.argv)
	a = ArgparseUi(parser)
	a.show()
	app.exec_()
	print ("Ok" if a.result() == 1 else "Cancel")
	if a.result() == 1: # Ok pressed
		parsed_args = a.parse_args()
		parse_files(parsed_args.swc, parsed_args.spine, parsed_args.output)
	else:
		args = None

def parse_files(swc_path, spine_path, output):
	dendrites = []
	with open(swc_path) as swc_file:
		reader = csv.reader(swc_file, delimiter=' ')

		cleaned_rows = []
		for row in reader:
			if not (len(row) <= 1 or row[0][0] == "#"):
				cleaned_rows.append(row)
		dendrites = list_dendrites_from_rows(cleaned_rows)

	with open(spine_path) as spine_file:
		reader = csv.reader(spine_file, skipinitialspace=True, delimiter=' ')
		first_row = True
		for row in reader:
			# Skip first row
			if first_row:
				first_row = False
			else:
				parent_id = int(row[13])
				for dendrite in dendrites:
					if dendrite.id == parent_id:
						dendrite.spine_count = dendrite.spine_count + 1

	valid_dendrites = filter_dendrites_based_on_soma_distance(dendrites, 60)
	remove_invalid_parents_and_children(valid_dendrites)
	# All dendrites contained in far_dendrite are further than 60 and their parent/children are contained in this array.
	branches = split2(dendrites, valid_dendrites)
	
	with open(output, 'w', newline='\n') as csvfile:
		writer = csv.writer(csvfile, delimiter=';')
		writer.writerow(["START","END","LENGTH","SPINE_COUNT","SPINE_DENSITY"])
		for branch in branches:
			total_spine_count = 0
			for dendrite in branch:
				total_spine_count = total_spine_count + dendrite.spine_count
			end_dendrite = find_branch_end(branch)
			start_dendrite = find_branch_start(branch)
			branch_length = calculate_branch_length(end_dendrite)

			writer.writerow([start_dendrite, end_dendrite, branch_length, total_spine_count, total_spine_count/branch_length])

def calculate_branch_length(branch_end):
	branch_length = 0.0
	current_dendrite = branch_end
	while(current_dendrite.parent is not None):
		branch_length = branch_length + current_dendrite.distance_with_parent()
		current_dendrite = current_dendrite.parent
	return branch_length

def find_branch_start(branch):
	for dendrite in branch:
		if dendrite.parent is None:
			return dendrite

def find_branch_end(branch):
	for dendrite in branch:
		if not dendrite.children:
			return dendrite

def remove_invalid_parents_and_children(dendrites):
	for dendrite in dendrites:
		if dendrite.parent not in dendrites:
			# This parent is not within the valid dendrite list
			dendrite.parent = None
		for children in dendrite.children:
			if children not in dendrites:
				# This children is not within the valid children list
				dendrite.children.remove(children)


def split2(dendrites, valid_dendrites):
	branches = []
	for dendrite in dendrites:
		if not dendrite.children:
			# We found an end!
			branch = visit_branch(dendrite, valid_dendrites)
			if branch:
				branches.append(branch)
	return branches

def visit_branch(end_dendrite, valid_dendrites):
	branch = []

	current_dendrite = end_dendrite
	while(current_dendrite in valid_dendrites):
		branch.append(current_dendrite)
		current_dendrite = current_dendrite.parent
	return branch

def split_into_branches(dendrites, branches=[]):

	if len(dendrites) > 0:
		seed_dendrite = dendrites[0]

		branch = [seed_dendrite]
		current_dendrite = seed_dendrite
		print("Starting children at {}".format(current_dendrite))
		while(len(current_dendrite.children) > 0):
			branch.append(current_dendrite.children[0]) # TODO - go through all children
			current_dendrite = current_dendrite.children[0]
		current_dendrite = seed_dendrite
		print("Starting parent at {}".format(current_dendrite))
		while(current_dendrite.parent is not None):
			branch.append(current_dendrite.parent)
			current_dendrite = current_dendrite.parent
		
		branches.append(branch)

		for dendrite in branch:
			print("Removing {}".format(dendrite))
			if dendrite in dendrites:
				dendrites.remove(dendrite)
		split_into_branches(dendrites, branches)

def filter_dendrites_based_on_soma_distance(dendrites, distance_in_micron):
	soma_count = 0
	soma = None
	for dendrite in dendrites:
		if dendrite.type == 1:
			soma_count = soma_count + 1
			soma = dendrite
	assert soma_count == 1

	valid_dendrites = []
	for dendrite in dendrites:
		if abs(dendrite.distance(soma)) >= distance_in_micron:
			valid_dendrites.append(dendrite)
	return valid_dendrites

def list_dendrites_from_rows(rows):
	dendrites = []
	for row in rows:

		parent_id = int(row[6])
		parent = None
		if not parent_id < 0:
			for dendrite in dendrites:
				if dendrite.id == parent_id:
					parent = dendrite

		dendrite = Dendrite(int(row[0]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), parent, int(row[1]))
		dendrites.append(dendrite)
	return dendrites

class Dendrite:
	def __init__(self, id,x,y,z,radius,parent=None, type=-1):
		self.id=id
		self.x=x
		self.y=y
		self.z=z
		self.radius=radius
		self.parent=parent
		self.children = []
		self.point = numpy.array((x,y,z))
		self.spine_count = 0
		self.type=type

		if parent is not None:
			parent.children.append(self)

	def distance(self, dendrite):
		return numpy.linalg.norm(self.point-dendrite.point)
	
	def distance_with_parent(self):
		if self.parent is None:
			return 0
		return self.distance(self.parent)

	def __str__(self):
		return "Dendrite#{}".format(self.id)

if __name__ == "__main__":
	main()