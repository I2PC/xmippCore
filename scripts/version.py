import argparse, os

def __getReleaseName(path: str, lineStart: str) -> str:
	"""
	### This function returns the line containing the given string in the provided file.

	#### Params:
	- path (str): Path to the file to be parsed.
	- lineStart (str): Start of the line to detect.

	#### Returns:
	- (str): First line found that starts with the given string.
	"""
	with open(path) as cmakeFile:
		for line in cmakeFile.readlines():
			if line.startswith(lineStart):
				return line[:-1] if line.endswith("\n") else line
	return ''

if __name__ == "__main__":
	# Generate and parse args
	parser = argparse.ArgumentParser(prog="version")
	parser.add_argument("-t", "--type", choices=['full', 'number', 'name'], default='full',
		help='Type of version to show.\n'
			'\'number\' only shows the version number.\n'
			'\'name\' only shows the release name.\n'
			'\'full\' shows the full release name, including number and name.')
	parser.add_argument("-k", "--keep-format", action="store_true")
	args = parser.parse_args()

	# Change directory to current file
	os.chdir(os.path.dirname(os.path.abspath(__file__)))

	# Get release full name
	relaseLine = __getReleaseName("../CHANGELOG.md", "##")
	fullName = relaseLine.replace("## Release ", "")

	# Print name deppending on argument provided
	resultName = fullName if args.keep_format else fullName.replace(" ", "")
	if args.type != 'full':
		parts = fullName.split("-")
		resultName = parts[0].strip() if args.type == 'number' else parts[1].strip()
	print(resultName)