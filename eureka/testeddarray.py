import ddarray

a = ddarray.double_ddarray(2, 3)
for i in range(2):
    for j in range(3):
        ddarray.setitem(a, i, j, i + 1)

print (ddarray.calculate(a, 2, 3))
for i in range(2):
	for j in range(3):
		print("{} ".format(ddarray.getitem(a,i,j)))
ddarray.deleteddarray(a,2,3)