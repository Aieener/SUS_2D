import numpy as np

data = np.loadtxt('dataplot.dat')
avetho = np.average(data[:10000, 4])
avever = np.average(data[:10000, 2])
avehor = np.average(data[:10000, 3])
aveQ = np.average(data[:10000, 1])
# read = "The average density Rho is ", avetho, '\n', "The average Vertical Rod is ", avever, \
#        '\n', "The average Horizontal Rod is ", avehor, '\n', "The average Up Rod is ", aveup
print "The average density Rho is %.10f \nThe average Vertical Rod is  %.5f" \
      "\nThe average Horizontal Rod is %.5f\nThe average Q is  %.5f" % (avetho,avever,avehor,aveQ)

with open("meanvalues.txt", "w") as file:
    file.write("The average density Rho is %.10f \nThe average Vertical Rod is  %.5f"
               "\nThe average Horizontal Rod is %.5f\nThe average Up Rod is  %.5f" % (avetho,avever,avehor,aveQ))
