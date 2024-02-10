import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random
import random
number=50
length=10

xyz_file=open('atoms_inbox.xyz','w')

coordinates={}


i=0
ditance=0
while i<number:
    x=random.random()*length
    y=random.random()*length
    z=random.random()*length
    coordinates[i]=[x,y,z]

    if i>0:
        distance_list=[]
        for j in range(0,i):
            dx=abs(coordinates[j][0]-coordinates[i][0])
            dy=abs(coordinates[j][1]-coordinates[i][1])
            dz=abs(coordinates[j][2]-coordinates[i][2])
#Using Minimum image convention
            if dx>length/2:
              dx=dx-length
            if dy>length/2:
                dy=dy-length
            if dz>length/2:
                dz=dz-length
#End of Minimum image convention
            distance=dx**2+dy**2+dz**2
            distance=distance**0.5
            distance_list.append(distance)

        if min(distance_list)>1:
            i+=1
    else:
        i+=1
xyz_file.write(f'{number}\n')
xyz_file.write(f'Atoms in box\n')
for atoms in range(0,number):
    xyz_file.write(f'H\t{format(coordinates[atoms][0], ".2f")}\t{format(coordinates[atoms][1], ".2f")}\t{format(coordinates[atoms][2], ".2f")}\n')
xyz_file.close()
print('XYZ file written')
#***********************************************************************************
# Read the coordinate file
def number_of_atoms(filename):
    coordinates_dict = {}
    x_coord=[]
    y_coord=[]
    z_coord=[]
    with open(filename, 'r') as file:
        num_atoms = int(file.readline().strip())
        comment = file.readline()
        atom_coordinates = file.readlines()
        i=0
        while i< num_atoms:
            values=atom_coordinates[i].split('\t')
            x_coord.append(float(values[1]))
            y_coord.append(float(values[2]))
            z_coord.append(float(values[3]))
            coordinates_dict[i]=[float(values[1]),float(values[2]),float(values[3])]
            i+=1
    return(int(num_atoms),coordinates_dict,x_coord,y_coord,z_coord)

# Calculate the energy of system
def calculate_energy(x_coord, y_coord, z_coord, length_of_box):
    length = length_of_box
    number_atoms = len(x_coord)
    tot_energy_initial=0
    # Create arrays for x, y, z coordinates
    x = np.array(x_coord)
    y = np.array(y_coord)
    z = np.array(z_coord)

    i=0
    while i<number_atoms:
        dx=abs(np.array(x[i:])-x[i])
        dx=dx[1:]
        dx[dx > length/2] = dx[dx > length/2]-length

        dy=abs(np.array(y[i:])-y[i])
        dy=dy[1:]
        dy[dy > length/2] = dy[dy > length/2]-length

        dz=abs(np.array(z[i:])-z[i])
        dz=dz[1:]
        dz[dz > length/2] = dz[dz > length/2]-length

        distance_squared = dx**2 + dy**2 + dz**2
        distance_6 = distance_squared ** 3
        distance_12 = distance_6 ** 2
        # Calculate energy
        tot_energy_initial += np.sum((1/distance_12 - 1/distance_6) * 4)
        i+=1
    return tot_energy_initial

#Montecarlo
def montecarlo(length,T,cycles,moves,delta,total_initial_energy,x_coord, y_coord, z_coord):
    i = 0
    potential_energy = []
    potential_energy.append(total_initial_energy)
    new_file = open('simulation.xyz', 'w')
    while i<cycles:
        j=0
        #print(f'Cycle: {i}')
        while j<moves:
            random_integer = random.randint(0, moves-1) #for picking up the particle
            values = [x_coord[random_integer],y_coord[random_integer],z_coord[random_integer]]
            x_coord[random_integer] = x_coord[random_integer]+(2*random.random()-1)*delta
            y_coord[random_integer] = y_coord[random_integer]+(2*random.random()-1)*delta
            z_coord[random_integer] = z_coord[random_integer]+(2*random.random()-1)*delta
            # Checking that atom should not go out of the box (0, length)
            if x_coord[random_integer] > length:
                x_coord[random_integer] = x_coord[random_integer] - length
            elif x_coord[random_integer] < 0:
                x_coord[random_integer] = x_coord[random_integer] + length
            if y_coord[random_integer] > length:
                y_coord[random_integer] = y_coord[random_integer] - length
            elif y_coord[random_integer] < 0:
                y_coord[random_integer] = y_coord[random_integer] + length
            if z_coord[random_integer] > length:
                z_coord[random_integer] = z_coord[random_integer] - length
            elif z_coord[random_integer] < 0:
                z_coord[random_integer] = z_coord[random_integer] + length
            new_energy=calculate_energy(x_coord, y_coord, z_coord, length)
            del_E=new_energy-total_initial_energy
            if del_E<0:
                total_initial_energy=new_energy
            else:
                calculation=np.exp(-del_E/T)
                metropolis1=random.random() #generating a random rumber
                #check if the random number is less than -del_E/T
                if calculation>metropolis1:
                    total_initial_energy=new_energy
                else:
                    x_coord[random_integer]=values[0]
                    y_coord[random_integer]=values[1]
                    z_coord[random_integer]=values[2]
            j+=1
        potential_energy.append(new_energy)
        if i%100==0:
            print(i)
            new_file.write(f'{number_atoms}\n')
            new_file.write(f'Atoms in box\n')
            for atoms in range(0, number_atoms):
                new_file.write(f'H\t{format(x_coord[atoms], ".2f")}\t{format(y_coord[atoms], ".2f")}\t{format(z_coord[atoms], ".2f")}\n')
        i+=1
    new_file.close()
    return(potential_energy)
print("DONE")
#Calling functions
length=10
T=0.50
cycles=100000
moves=50
delta=0.05 #move the particle by this value

start_time = time.time()
number_atoms, coordinates,x_coord,y_coord,z_coord=number_of_atoms('atoms_inbox.xyz')
print(f'Number of atoms:{number_atoms}')
total_initial_energy = calculate_energy(x_coord, y_coord, z_coord, length)
potential_energy=montecarlo(length,T,cycles,moves,delta,total_initial_energy,x_coord, y_coord, z_coord)
end_time = time.time()

plt.plot(np.arange(1,len(potential_energy)+1),potential_energy)
plt.xlabel('Number of steps')
plt.ylabel('Total energy ')
plt.savefig("potential_energy.png", dpi=100)
print(f'Time taken: {end_time - start_time} seconds')
