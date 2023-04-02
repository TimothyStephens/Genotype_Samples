#!/usr/bin/env python
DESCRIPTION = '''
add_value_to_table - Will add values into a table (via a new end column)
                     given key:value pairs. 

NOTE:
	- Ignore comment ('#' by defult; can be turned off) and blank lines	
	- Uses exact string matching. Can not do regex or partial matching
	- Assumes first column is 'key' and collowing column/s are 'value'
	- Assumes key is unique in -a/--add; if not unique will take the last value and print warning.
	- Will always add 'blank' values if target column not in key:value pairs.

Uses SQLite3 to work quickly with large --add files.

'''
import sys
import argparse
import logging
import sqlite3

VERSION=0.1

## Pass arguments.
def main():
	# Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-i', '--input', metavar='data_file.txt', default=sys.stdin, type=argparse.FileType('r'), required=False, help='Input file (default: stdin)')
	parser.add_argument('-o', '--output', metavar='data_file_with_extra_column.txt', default=sys.stdout, type=argparse.FileType('w'), required=False, help='Output file (default: stdout)')
	parser.add_argument('-a', '--add', metavar='info_to_add.txt', type=argparse.FileType('r'), required=True, help='Input key:value pairs')
	parser.add_argument('-c', '--col', default=1, type=int, required=False, help='Column in --input of interest')
	parser.add_argument('-d', '--default', default='', type=str, required=False, help='Value to add if not in -a/--add (default: %(default)s)')
	parser.add_argument('--delim_input', default='\t', type=str, required=False, help='Delimiter for --input (default: \\t)')
	parser.add_argument('--delim_add', default='\t', type=str, required=False, help='Delimiter for --add (default: \\t)')
	parser.add_argument('--keep_comments', action='store_true', required=False, help='Keep comment lines from input file in output file (default: %(default)s)')
	parser.add_argument('--debug', action='store_true', required=False, help='Print DEBUG info (default: %(default)s)')
	args = parser.parse_args()
	
	# Set up basic debugger
	logFormat = "[%(levelname)s]: %(message)s"
	logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	
	logging.debug('%s', args) ## DEBUG
	
	db = sqlite3.connect(":memory:")
	load_key_value_from_file(db, args.add, args.delim_add)
	add_new_column(db, args.input, args.output, args.col, args.default, args.delim_input, args.keep_comments)



def add_new_column(db, input_file, output_file, col, default, delim_input, keep_comments):
	c = db.cursor()
	# For each line in input file
	for line in input_file:
		line = line.strip('\n')
		if not line:
			continue
		
		if line.startswith('#'):
			if keep_comments:
				output_file.write(line + '\n')
			continue
		
		line_sep = line.split(delim_input)
		try:
			key = line_sep[col-1]
		except IndexError:
			logging.info("[ERROR]: %s", line)
			logging.info("[ERROR]: -c/--col %s out of range for --infile", col)
			sys.exit(1)
		
		c.execute('SELECT value FROM values2add WHERE key=?', (key, ))
		value = [x[0] for x in c.fetchall()]
		logging.debug('Values selected for key %s : %s', key, value) ## DEBUG
		
		if len(value) > 1:
			logging.warning('%s occurs multiple times - taking just one entry.', key) ## DEBUG
			output_file.write(line + delim_input + value[0] + '\n')
			logging.debug('Value added: %s:%s', key, value[0]) ## DEBUG
		elif len(value) == 1:
			output_file.write(line + delim_input + value[0] + '\n')
			logging.debug('Value added: %s:%s', key, value[0]) ## DEBUG
		else:
			output_file.write(line + delim_input + default + '\n')
			logging.debug('Default added: %s', default) ## DEBUG



def load_key_value_from_file(db, keyvalue_file, delim):
	'''
	Loads a table of key:value pairs using SQLite3.
	'''
	c = db.cursor()
	c.execute('CREATE TABLE values2add (key, value)')
	c.execute('CREATE INDEX key_index ON values2add (key)')
	db.commit()
	
	c = db.cursor()
	for line in keyvalue_file:
		line = line.strip()
		if not line or line.startswith('#'):
			continue
		
		line_split = line.split(delim)
		key = line_split[0]
		value = delim.join(line_split[1:])
		c.execute('INSERT INTO values2add (key, value) VALUES (?, ?)', (key, value))
	db.commit()


if __name__ == '__main__':
	main()
