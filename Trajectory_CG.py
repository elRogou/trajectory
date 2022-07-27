#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import re


# In[2]:


class atom:
    """
    Atributes::  index, name (CA in coarse grained simulations),
    and aa it represents.
    ---------------------------------
    """
    __chrg_lookup = {'ARG': +1,
                   'LYS': +1,
                   'ASP': -1,
                   'GLU': -1,
                   'HIS': 0,
                   'SER': 0,
                   'THR': 0,
                   'ASN': 0,
                   'GLN': 0,
                   'CYS': 0,
                   'GLY': 0,
                   'PRO': 0,
                   'ALA': 0,
                   'ILE': 0,
                   'LEU': 0,
                   'MET': 0,
                   'PHE': 0,
                   'TRP': 0,
                   'TYR': 0,
                   'VAL': 0,
                   'P'  : 0,
                   'B'  : 0,
                   'S'  : 0,
                   'DG' : 0,
                   'DC' : 0,
                   'DT' : 0,
                   'DA' : 0,
                   'G'  : 0,
                   'C'  : 0,
                   'T'  : 0,
                   'A'  : 0
                  }
    
    def __init__(self,index, name, amino_acid):
        self.index     = int(index)
        self.name      = name
        self.aa        = amino_acid
        self.charge    = self.__chrg_lookup[self.aa.upper()]
        
    def __str__(self):
        return f'{self.index} {self.name} {self.aa}\n'


# In[3]:


class structure:
    """
    Atributes:: n_atoms, atoms (dictionary)
    ---------------------------------
    atoms can be added to the dictionary by the append method
    """
    def __init__(self,n_atoms):
        self.n_atoms = n_atoms
        self.atoms   = {}
        
    def append(self,atom, chain):
        if chain not in self.atoms.keys():
            self.atoms[chain] = {}
        self.atoms[chain][atom.index] = atom
        
    def totalCharge(self):
        totalQ = 0
        for chain in self.atoms.values():
            for atom in chain.values():
                totalQ += atom.charge
        return totalQ


# In[1]:


class coordinate:
    """
    Atributes:: x, y, z, pos
    ---------------------------------
    method: measuring distance between two coordinates
    """
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
        self.pos = np.array([x,y,z])
        
    def __str__(self):
        return f'{self.x}\t{self.y}\t{self.z}\n'
    
    def get_distance(self, coord):
        distance = np.linalg.norm(self.pos - coord.pos) #in nm
        return distance


# In[5]:


class step:
    """
    Atributes:: n_step, coordinates
    ---------------------------------
    holds numpy array with the coordinates
    for all atoms in the sumulation
    """
    def __init__(self,n_step, n_atoms):
        self.n = n_step
        self.coordinates = np.zeros(n_atoms, dtype = object)
        
    def append(self, coord, i):
        self.coordinates[i] = coord


# In[6]:


class timesteps:
    """
    Atributes:: n_steps, steps
    ---------------------------------
    holds the total number of steps, numpy array with the 
    objects of class step
    """
    step_jump = 1000
    def __init__(self,n_steps):
        self.n_steps = n_steps
        self.steps = np.zeros(n_steps // self.step_jump, dtype = object)
    
    def append(self,step):
        self.steps[(step.n // self.step_jump) - 1] = step


# In[7]:


class Trajectory:
    """
    Atributes:: structure, timesteps
    ---------------------------------
    The trajectory class object holds structure class (basically atoms)
    and the timesteps with coordinates at each step for each atom
    """
    def __init__(self,n_atoms, n_steps, temperature):
        self.structure   = structure(n_atoms)
        self.timesteps   = timesteps(n_steps)
        self.temperature = temperature
        
    def print(self):
        print(f'The following trajectory contains {self.structure.n_atoms} atoms\n            The number of steps is {self.timesteps.n_steps/self.timesteps.step_jump:.1E}\n            Total simulation time is {self.timesteps.n_steps:.1E}            The temperature in this run was {self.temperature}\n'
            )


# In[10]:


class TrajectoryFile:
    """
    Attributes:: filename
    ---------------------------------
    public methods:: createTrajectoryObject reads the file
    and creates object of class Trajectory
    """
    __steps_re = 'end\n(.*)\n'
    
    def __init__(self,filename):
        self.filename = filename
        
    def createTrajectoryObject(self):
        print(' Creating trajectory '.center(40,'='))
        with open(self.filename) as trj:
            n_atoms, n_steps, n_chain = self.__getAtomsAttributes(trj)
            temperature = self.__getTemperature()
            print(f'The number of atoms is {n_atoms} in one chain\n')
            print(f'The number of steps in the trajectory is {n_steps}\n')
            print(f'The temperature in the trajectory is {temperature}\n')
            self.__Traj = Trajectory(n_atoms, n_steps, temperature)
            for i in range(n_chain+1):
                line = trj.readline()
            self.__getAtoms(trj, n_atoms, n_chain)
            line = trj.readline()
            assert str(n_atoms) in line
            self.__getSteps(trj,n_steps,n_atoms)
        return self.__Traj
    
    def __getTemperature(self):
        split = self.filename.split('_')
        temperature = float(split[-2])
        return temperature
        
    def __getAtoms(self, trj, n_atoms, n_chain):
        print(n_atoms, n_chain)
        for i in range(n_atoms):
            line = trj.readline()
            split_line = line.split()
#             print(split_line)
            a = atom(*split_line)
            self.__Traj.structure.append(a,n_chain)
                
    def __getSteps(self, trj, n_steps,n_atoms):
        line = trj.readline()
        n = int(line)
        current_step = step(n, n_atoms)

        #getCoords
        for i in range(n_atoms):
            line = trj.readline()
            split_line = re.findall('.{1,8}', line)
            try:
                coord = coordinate(*map(float,split_line))
            except:
                coord = coordinate(0.0,0.0,0.0)
                print(f'did not handle {i+1} atom coord')
                
            current_step.append(coord,i)
        self.__Traj.timesteps.append(current_step)
#         print(current_step.n, n_steps)
        
        if current_step.n == n_steps:
            print(' Returning trajectory '.center(50,'='))
        else:
            trj.readline()
            Traj = self.__getSteps(trj, n_steps, n_atoms)
            
    def __getAtomsAttributes(self, trj):
        n_atoms = 0
        n_chain = int(trj.readline())
        for i in range(n_chain):
            n_atoms += int(trj.readline())
        n_steps = self.__findNumberOfSteps(trj)
        return  n_atoms, n_steps, n_chain
    
    def __findNumberOfSteps(self, trj):
        match = re.search(self.__steps_re, trj.read())
        if match:
            n_steps = int(match.group(1))
            trj.seek(0)
            return n_steps
        else:
            print('did not find the number of steps')


# In[ ]:




