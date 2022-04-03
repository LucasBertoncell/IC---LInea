import numpy as np
import matplotlib.pyplot as plt
from sys import *
from os import *


################################## Percentual Difference #####################################

## getting values

print "\nTwo files containing xi (correlation values) will be compared\n"

user_input_1 = raw_input("Enter the path to the file 1: ")
assert os.path.exists(user_input_1), " \nFile not found at: "+str(user_input)
f = open(user_input_1,'r')  #r stands for reading

xi_1 = user_input_1[:, 1]  # second column stored in array 1
separation_1 = user_input_1[:, 0] # can be theta, r, or any type of CF

f.close()


user_input_2 = raw_input("\nEnter the path to the file 2: ")
assert os.path.exists(user_input_2), " \nFile not found at: " +str(user_input)

f = open(user_input_2,'r')  #r stands for reading


xi_2 = user_input_1[:, 1]  # second column stored in array 1
separation_2 = user_input_1[:, 0]

f.close()


## getting percentual difference
diff = []
i = 0
while i < len(xi_1):

	diff = ((xi_2  - xi_1)/xi_1)*100 #array
	i += 1
## plotting

option = raw_input("\nPlot percentual diference vs separation? (y/n):")

if option == y:

	plt.figure(figsize = (7,6), dpi = 2048)
	plt.subplot(311)
	plt.plot(separation, diff, 'r-', label="Percentual Difference")
	plt.xlabel("$\Theta$ [deg]")
	plt.ylabel("percentual difference [%]")

	plt.subplot(312)
	plt.plot(separation_1, xi_1, 'g', label="xi_1")
	plt.xlabel("$\Theta$ [deg]")
	plt.ylabel("$w(\Theta)$")

	plt.subplot(313)
	plt.plot(separation_2, xi_2, 'b', label="xi_2")
	plt.xlabel("$\Theta$ [deg]")
	plt.ylabel("$w(\Theta)$")
	plt.savefig('Percentual_Difference.png')
	plt.show()

else:

	sys.exit()