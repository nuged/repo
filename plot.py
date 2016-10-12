import matplotlib.pyplot as plt

gd = open("gd.txt","r")
for line in gd.readlines():
    y = line.split()
    x = [a for a in range(len(y))]
    plt.plot(x, y)
plt.xlim(0, 1000)
plt.show()
gd.close()

met = open("met.txt", "r")
for line in met.readlines():
    y = line.split()
    x = [a for a in range(len(y))]
    plt.plot(x, y)
plt.xlim(0, 1000)
plt.show()
met.close()

ann = open("ann.txt", "r")
for line in ann.readlines():
    y = line.split()
    x = [a for a in range(len(y))]
    plt.plot(x, y)
plt.xlim(0, 1000)
plt.show()
ann.close()

