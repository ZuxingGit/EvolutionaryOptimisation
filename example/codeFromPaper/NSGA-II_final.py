# https://www.sciencedirect.com/science/article/pii/S2352710224004108#cebib0010
# Multi-objective optimization of a residential zone by proposing appropriate comfort factors using none dominated sorting genetic algorithm
#Libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from copy import deepcopy
from random import sample
from time import sleep
#constants
'''
gene1 = width
gene2 = height
gene3 = rotation
gene4 = WWR
gene5 = V-shader01
gene6 = V-shader02
gene7 = H-shader
'''
num_genes = 7
domain_genes = [0,1]
num_obj = 3

#Hyperparameters
num_population = 34
num_cross = 24
num_mutation = 8
num_generation = 20
alpha = 0.5
beta = 10




#class definition
class Chromosom:
    population = []
    
    def __init__(self,genome=None):
        if genome:
            self.genome = np.array(genome)
        else:
            self.genome = np.random.uniform(domain_genes[0],domain_genes[1],num_genes)
        self.fitness = self.bridge()
        Chromosom.population.append(self)
    
    def bridge(self):
        df = pd.DataFrame(self.genome)
        df.to_csv('for_grass.csv')
        print('write')
        sleep(300)
        print('read')
        f = pd.read_csv('for_python.csv')
        print(f)
        print(np.array([f.values[0][0],f.values[0][1],f.values[0][2]]))
        return np.array([f.values[0][0],f.values[0][1],f.values[0][2]])
    
    def blend_crossover(self,other):
        landa = np.random.uniform(-alpha, 1+alpha, num_genes)
        G1 = self.genome
        G2 = other.genome
        offspring1 = landa*G1 + (1-landa)*G2
        offspring2 = landa*G2 + (1-landa)*G1
        
        offspring1[offspring1>domain_genes[1]] = domain_genes[1]
        offspring1[offspring1<domain_genes[0]] = domain_genes[0]
        offspring2[offspring2>domain_genes[1]] = domain_genes[1]
        offspring2[offspring2<domain_genes[0]] = domain_genes[0]
        
        [Chromosom(list(offspring1)), Chromosom(list(offspring2))]
    
    
    def mutate(self):
        temp_genome = deepcopy(self.genome)
        sigma = (domain_genes[1] - domain_genes[0]) / beta
        temp_genome[np.random.randint(0,num_genes)] += np.random.normal(0,sigma)
        
        temp_genome[temp_genome>domain_genes[1]] = domain_genes[1]
        temp_genome[temp_genome<domain_genes[0]] = domain_genes[0]
        
        Chromosom(list(temp_genome))
    
    
    @classmethod
    def non_dominated_sorting(cls):
        counter = 1
        n_d_s = {'front1':[]}
        num = len(cls.population)
        
        for i in cls.population:
            i.dominated_count = 0
            i.dominated_chromes = []
            i.grade = None
        
        for i,j in combinations(range(num),2):
            chrome1 = cls.population[i]
            chrome2 = cls.population[j]
            if all(chrome1.fitness<=chrome2.fitness) and (any(chrome1.fitness<chrome2.fitness)):
                chrome1.dominated_chromes.append(j)
                chrome2.dominated_count += 1
            
            elif all(chrome2.fitness<=chrome1.fitness) and (any(chrome2.fitness<chrome1.fitness)):
                chrome2.dominated_chromes.append(i)
                chrome1.dominated_count += 1
        #for i in range(num):
            #print(cls.population[i].dominated_count)
        
        for i in range(num):
            if cls.population[i].dominated_count == 0:
                n_d_s['front1'].append(i)
                cls.population[i].grade = 1
                #print(i)
                #print(cls.population[i].grade)
        #print(n_d_s['front1'])
        #print('-----------------------')
        while True:
            counter += 1
            n_d_s[f'front{counter}'] = []
            
            for i in n_d_s[f'front{counter-1}']:
                chrome1 = cls.population[i]
                for j in chrome1.dominated_chromes:
                    chrome2 = cls.population[j]
                    chrome2.dominated_count -= 1
                    if chrome2.dominated_count == 0:
                        n_d_s[f'front{counter}'].append(j)
                        chrome2.grade = counter
            if len(n_d_s[f'front{counter}']) == 0:
                del n_d_s[f'front{counter}']
                break
        #print(n_d_s)
        return n_d_s
    
    
    @classmethod        
    def crowding(cls,fronts):
        
        for i in cls.population:
            i.crowding_distance = 0
            
        for front in fronts.values():
            n = len(front)
            if n>1:
                costs = np.vstack([cls.population[t].fitness for t in front])
                d = np.zeros((n,num_obj))
                
                for j in range(num_obj):
                    c = np.sort(costs[:,j])
                    print(c)
                    idx = np.argsort(costs[:,j])
                    
                    d[idx[0]][j] = np.inf
                    d[idx[-1]][j] = np.inf
                    
                    for i in range(1,n-1):
                        d[idx[i]][j] = np.abs(c[i+1]-c[i-1]) / np.abs(c[-1] - c[0])
                
                for i,j in enumerate(front):
                    cls.population[j].crowding_distance = np.sum(d[i,:])
            
            else:
                cls.population[front[0]].crowding_distance = np.inf
                
    
    @classmethod
    def final_sort(cls):
        cls.population = sorted(cls.population, key=lambda t: t.crowding_distance, reverse=True)
        cls.population = sorted(cls.population, key=lambda t: t.grade)      
        

    @classmethod
    def create(cls):
        for _ in range(num_population*4):
            cls()
            
        for it in range(num_generation):
            print(it)
            mom = sample(cls.population,num_cross)
            for i in range(0,num_cross,2):
                parent1 = mom[i]
                parent2 = mom[i+1]
                parent1.blend_crossover(parent2)
            
            
            son = sample(cls.population,num_mutation)
            for i in son:
                i.mutate()
                
            n_d_f = cls.non_dominated_sorting()
            cls.crowding(n_d_f)
            cls.final_sort()
            cls.population = cls.population[:num_population]
            
                
            ff = np.array([i.fitness for i in cls.population if i.grade == 1])
            ff = np.array(sorted(ff,key=lambda t: t[0]))
            

            
            df = {'f1': ff[:,0],'f2': ff[:,1],'f3': ff[:,2]}
            final_res = pd.DataFrame(df)
            final_res.to_csv('firstfront.csv')
            print(final_res)
Chromosom.create()
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            