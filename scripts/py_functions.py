#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-14
# Last Modified: 2020-04-14



"""  """

# timer
import time
start = time.time() # start the timer from import

def timer():
	
	end = time.time()
	duration = end-start
	duration = round(duration)
	string = "\n..........................\n"\
			"   Time elapsed: {} sec"\
			"\n..........................\n"\
			.format(duration)

	print(string)
