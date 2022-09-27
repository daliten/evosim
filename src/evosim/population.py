import numpy as np

def rescale_mut(tstep,mut_num):
    rescaled_mut = mut_num*tstep
    if int(rescaled_mut)<1:
        num_start = np.random.binomial(mut_num,rescaled_mut)
        if num_start>0:
            return_dict={'start':num_start,'end':0}
        else:
            return_dict={'start':0,'end':mut_num}
    else:
        return_dict={'start':int(rescaled_mut),'end':0}
    return return_dict

class Population:
     """
     Defines a Population object and provides functions for simulating a continuous-time branching process through either extended-time leaping or Gillespie's algorithm.
     
     """
     def __init__(self):
          
        self.pheno_dict = {}
        self.mut_dict = {}
        self.birth_rates = {}
        self.death_rates = {}
        self.capacity = 0.
        
        self.leap_thresh = 10**2
                
        self.t=0


     def reset(self,pop_dict,**kwargs):
         """
         

         Parameters
         ----------
         pop_dict : dict
             Initial population distribution specified as {'Name1':Pop1,'Name2':Pop2,...,'NameN':PopN} with Pop1,...,PopN int; 'Name1',...,'NameN' can be any user-specified type.
         leap_threshold: int, optional 
             Subpopulation threshold at which leaping sets in, if in 'leaping' mode. Defaults to 10^2.
         carrying_capacity: int, optional 
             Carrying capacity of environment. Defaults to 10^18.

         Returns
         -------
         None.

         """

         self.pheno_dict = pop_dict #for cell types 1,2,3 this is {1:pop1,2:pop2,3:pop3}
         
         self.leap_thresh = kwargs.get('leap_threshold',10**2)
         self.capacity = min(kwargs.get('carrying_capacity',10**18),10**18)

         self.t=0
     


     
     def get_char_time(self,key):
         """
         Computes the current characteristic times for a single cell type.

         Parameters
         ----------
         key : user-specified type
             Cell type identifier (as defined in self.pheno_dict).
         pop: int
             Current population of this cell type.
             
         Raises
         ------
         Exception
            If there is a species for which neither birth nor death rate is positive.
            
         Returns
         -------
         float
             Characteristic time.

         """

         
         [br,dr]=self.birth_death_rates(key)
         pop=self.pheno_dict[key]
         max_rate = max(br,dr)
         if max_rate > 0:
             t_char = 1/(pop*max_rate)
         else:
             raise Exception('At least one of the birth or death rate must be positive for each cell type.')
         if pop<10:
             t_char = t_char/2
         return t_char

     
     def get_pops(self):
        """

         Returns
         -------
         list
             List of current population levels, ordered in accordance with types specified in self.pheno_dict.

         """
        return list(self.pheno_dict.values())

    
     def birth_death_rates(self,key):
         """
         Computes the current per-cell birth and death rates for a specified cell type. Assumes: logistic growth; resource capacity only affects birth rates.

         Parameters
         ----------
         key : user-specified type
             Cell type identifier (as defined in self.pheno_dict).

         Returns
         -------
         list
             A two-element list [birth rate, death rate].

         """
         b_rate = self.birth_rates[key]*np.clip((1-sum(self.get_pops())/self.capacity),0,1) #birth rate for logistic growth
         d_rate = self.death_rates[key]  
         return [b_rate,d_rate]

     def r_rates(self,key):
         """
         Computes the per-population birth and death rates for a specified cell type .

         Parameters
         ----------
         key : user-specified type
             Cell type identifier (as defined in self.pheno_dict).

         Returns
         -------
         list
             A two-element list [birth rate, death rate].

         """
         this_pop = self.pheno_dict[key]
         [b_rate,d_rate]=self.birth_death_rates(key)
         return [this_pop*b_rate,this_pop*d_rate]

             
     def eff_branch_times(self):
         """
         Raises
         ------
         Exception
            If both the birth and death rates for one of the cell types are zero.
            
         Returns
         -------
         float
             Two dictionaries: one of natural branching times (pops below leap threshold) and one of effective branching times (pops above leap threshold)

         """
         branch_times={}
         for key, pop in self.pheno_dict.items():
             if pop>0:
                 [br,dr]=self.birth_death_rates(key)
                 max_rate = max(br,dr)
                 if max_rate == 0:
                     raise Exception('At least one of the birth or death rate must be positive for each cell type.')
                 branch_times[key]=1/(min(pop,self.leap_thresh)*max_rate)
             
         return branch_times


             
     def branching_sim_multinomial(self,key,tstep):
         """
         Performs a single extended-time leaping simulation for a given cell type over a specified time interval.
         First computes the rescaled probabilities for this population level and time interval; then uses probabilities to drawn from a multinomial distribution.
         Updates self.pheno_dict with new population structure.

         Parameters
         ----------
         key : user-specified type
             Cell type identifier (as defined in self.pheno_dict).
         tstep : float
             Time interval.

         Returns
         -------
         None.

         """
         
         pop_level = self.pheno_dict[key]
         [br,dr]=self.birth_death_rates(key)
         

         tchar=self.get_char_time(key) #compute the characteristic time for this population

         if tchar<tstep: #rescale br,dr if condition holds         
            if br>=dr:
                br_base = br
                br=(np.exp(tstep*(br-dr))-1)/tstep+dr
                prob_mut_factor = (np.exp(tstep*(br_base-dr))-1)/((1-dr/br_base)*(np.exp(tstep*(br_base-dr))-1+dr*tstep))
            else:
                dr_base = dr
                dr=(1-np.exp(tstep*(br-dr)))/tstep+br
                prob_mut_factor=(np.exp(tstep*(br-dr_base))-1)/((br-dr_base)*tstep)


         #Compute the probability vector for multinomial sampling
         if len(self.mut_dict)>0:
             mut_keys = self.mut_dict[key]
             mut_probs = list(self.mut_dict[key].values())
             if tchar<tstep:
                 mut_probs = np.array(mut_probs)*prob_mut_factor
             br_no_mut_prob = br*(1-sum(mut_probs))*tstep
             dr_prob = dr*tstep
             individual_br_mut_probs = [br*mut_probs[i]*tstep for i in range(len(mut_probs))]         
             no_change_prob = 1-(br_no_mut_prob+sum(individual_br_mut_probs)+dr_prob)
             self_probs = [no_change_prob, br_no_mut_prob,dr_prob]         
             prob_array=self_probs+individual_br_mut_probs
         else:
             br_no_mut_prob = br*tstep
             dr_prob = dr*tstep       
             no_change_prob = 1-(br_no_mut_prob+dr_prob)
             self_probs = [no_change_prob, br_no_mut_prob,dr_prob]
             prob_array=self_probs
         if any(p < 0 for p in prob_array):
             print(prob_array)
             print('tstep = ',tstep)
             print('tchar = ',tchar)
             print('pop = ',pop_level)
             print('br = ',br)
             print('dr = ',dr)
             print('br_no_mut_prob = ',br_no_mut_prob)


         #Multinomial sampling
         status = np.random.multinomial(pop_level,prob_array)

         #The new population level is the sum of cells sampled to undergo no change and twice the number of cells sampled to divide
         new_pop_level=max(0,status[0]+2*status[1])
     
         #Update the population level
         self.pheno_dict[key] = new_pop_level

         #Update the population levels from any mutations
         if len(self.mut_dict)>0:
             mut_status = status[3::]
             mut_keys = list(self.mut_dict[key].keys())
             for mut_ind in range(len(mut_status)):             
                 mut_key = mut_keys[mut_ind]
                 if mut_status[mut_ind]>0:   
                     self.pheno_dict[mut_key] = self.pheno_dict[mut_key] + mut_status[mut_ind]
               
                     
     def evolve_leaping(self,overall_time,sim_time,birth_rates,death_rates,mut_probs,**kwargs):
         """
         Evolves the system under specified parameters with the extended-time leaping algoritm.

         Parameters
         ----------
         overall_time : float
             Current time; time for starting simulation.
         sim_time : float
             Simulation time.
         birth_rates : dict
             Birth/division rates (inverse expected time to the division of a cell) specified as {'Name1':birth_rate1,'Name2':birth_rate2,...,'NameN':birth_rateN} with birth_rate1,...,birth_rateN floats; 'Name1',...,'NameN' can be any user-specified type.
         death_rates : dict
             Death rates (inverse expected time to the death of a cell) specified as {'Name1':death_rate1,'Name2':death_rate2,...,'NameN':death_rateN} with death_rate1,...,death_rateN floats; 'Name1',...,'NameN' can be any user-specified type.
         mut_probs : dict
             Mutation probabilities (probability of mutation to phenotype j in a single division of a phenotype i cell) specified as a dictionary of embedded dictionaries {'Name1':{'Name2':prob_mut_1_to_2,...,'NameN':prob_mut_1_to_N},...,'NameN':{'Name1':prob_mut_N_to_1,...,'Name(N-1)':prob_mut_N_to_(N-1)}}.
         detection_threshold: int, optional
             Stops the simulation when the population at or above specified level is detected. Default is no detection threshold (-1). Note that since population size is only checked at certain intervals during the simulation, threshold detection is only approximate.
         quit_at_zero_cells: bool, optional
             Stops the simulation when the total population size is zero. Default is True.
            

         Returns
         -------
         float
             Current updated time.
         is_done : bool
             Indicator for whether the simulation has terminated.

         """

         is_done = False         
         detection_pop = kwargs.get('detection_threshold',-1)
         quit_at_zero_cells=kwargs.get('quit_at_zero_cells',True)

         self.birth_rates = birth_rates 
         self.death_rates = death_rates
         self.mut_dict = mut_probs

         if sum(self.get_pops())<1: #no cells left, must quit sim
             if quit_at_zero_cells:
                 is_done = True
                 return overall_time,is_done
             else:
                 is_done = False
                 return overall_time+sim_time,is_done
         
            
         self.t = 0         
         while (self.t<sim_time and sum(self.get_pops())>0): #loop until either simulation time achieved or population drops to zero

             if sum(self.get_pops())<1:
                 break 
             
             while self.t<sim_time and sum(self.get_pops())>0:
                 #compute all effective branching times      
                 branch_times = self.eff_branch_times()
                 ordered_branch_times = dict(sorted(branch_times.items(), key=lambda item: item[1]))
                 #find the minimal branching time
                 tstep_branching = min(list(branch_times.values()))
                 if sum(self.get_pops())<10:
                     tstep_branching = tstep_branching/2
                 tstep_branching=min(tstep_branching,sim_time-self.t) #make sure not to exceed the actual remaining simulation time
                 
                 for k in ordered_branch_times.keys():  #k is pop key, going from small to large branching times
                     self.branching_sim_multinomial(k,tstep_branching)
                 self.t+=tstep_branching
                 if detection_pop>0 and sum(self.get_pops())>=detection_pop: #applies only to sims for population growth where a detection limit was specified. Quitting if detection limit achieved. 
                     is_done = True
                     return overall_time+self.t,is_done
                
             for k,v in self.pheno_dict.items(): 
                 if v<0:
                     self.pheno_dict[k]=0

             if sum(self.get_pops())<1: #no cells left, must quit sim
                 if quit_at_zero_cells:
                     is_done = True
                     return overall_time,is_done
                 else:
                     is_done = False
                     return overall_time+sim_time,is_done
             
             if detection_pop>0 and sum(self.get_pops())>=detection_pop: #applies only to sims for population growth where a detection limit was specified. Quitting if detection limit achieved. 
                 is_done = True
                 return overall_time+self.t,is_done
        
         return overall_time+self.t,is_done           
                
         
     def evolve_SSA(self,overall_time,sim_time,birth_rates,death_rates,mut_probs,**kwargs): 
         """
         Evolves the system under specified parameters with the SSA / Gillespie's algoritm.

         Parameters
         ----------
         overall_time : float
             Current time; time for starting simulation.
         sim_time : float
             Simulation time.
         birth_rates : dict
             Birth/division rates (inverse expected time to the division of a cell) specified as {'Name1':birth_rate1,'Name2':birth_rate2,...,'NameN':birth_rateN} with birth_rate1,...,birth_rateN floats; 'Name1',...,'NameN' can be any user-specified type.
         death_rates : dict
             Death rates (inverse expected time to the death of a cell) specified as {'Name1':death_rate1,'Name2':death_rate2,...,'NameN':death_rateN} with death_rate1,...,death_rateN floats; 'Name1',...,'NameN' can be any user-specified type.
         mut_probs : dict
             Mutation probabilities (probability of mutation to phenotype j in a single division of a phenotype i cell) specified as a dictionary of embedded dictionaries {'Name1':{'Name2':prob_mut_1_to_2,...,'NameN':prob_mut_1_to_N},...,'NameN':{'Name1':prob_mut_N_to_1,...,'Name(N-1)':prob_mut_N_to_(N-1)}}.
         detection_threshold: int, optional
             Stops the simulation when the population at or above specified level is detected. Default is no detection threshold (-1). Note that since population size is only checked at certain intervals during the simulation, threshold detection is only approximate.
         quit_at_zero_cells: bool, optional
             Stops the simulation when the total population size is zero. Default is True.
            

         Returns
         -------
         float
             Current updated time.
         is_done : bool
             Indicator for whether the simulation has terminated.


         """

         is_done = False
         detection_pop = kwargs.get('detection_threshold',-1)
         quit_at_zero_cells=kwargs.get('quit_at_zero_cells',True)
         
         self.birth_rates = birth_rates 
         self.death_rates = death_rates
         self.mut_dict = mut_probs

         if sum(self.get_pops())<1:
             if quit_at_zero_cells:
                 is_done = True
                 return overall_time,is_done
             else:
                 is_done = False
                 return overall_time+sim_time,is_done         

         tot_t = 0   
         while (tot_t<sim_time and sum(self.get_pops())>0):
             #self.step_counter += 1
             all_r = {}
             for key, pop in self.pheno_dict.items():
                 if pop>0:
                     [b,d]=self.r_rates(key)
                     if len(self.mut_dict)>0:
                         mut_probs_dict=self.mut_dict[key]
                         all_r[(key,'b')]=b*(1-sum(list(mut_probs_dict.values())))
                         all_r[(key,'d')]=d                 
                         for mut_key in mut_probs_dict.keys():
                             all_r[(key,mut_key)]=b*mut_probs_dict[mut_key]
                     else:
                         all_r[(key,'b')]=b
                         all_r[(key,'d')]=d           

             r_sum = sum(list(all_r.values()))
             all_r = {k: v/r_sum for k, v in sorted(all_r.items(), key=lambda item: item[1])}

             [r1,r2]=np.random.uniform(0,1,2)
             t = np.log(1/r1)/r_sum
             cumul_sum=0
             rxn_key='none'
             
             for k,v in all_r.items():
                 cumul_sum = all_r[k]+cumul_sum
                 if r2<=cumul_sum:
                     rxn_key = k
                     break

             if rxn_key!='none':
                 if rxn_key[1] == 'b': #division
                     self.pheno_dict[rxn_key[0]]=self.pheno_dict[rxn_key[0]]+1
                 elif rxn_key[1] == 'd': #death
                     self.pheno_dict[rxn_key[0]]=self.pheno_dict[rxn_key[0]]-1
                 else: #mutation
                     self.pheno_dict[rxn_key[1]]=self.pheno_dict[rxn_key[1]]+1

             for k,v in self.pheno_dict.items():
                 if v<0:
                     self.pheno_dict[k]=0
                                         
             tot_t += t
             
             if tot_t == 0:
                 break

             if sum(self.get_pops())<1:
                 if quit_at_zero_cells:
                     is_done = True
                     return overall_time,is_done
                 else:
                     is_done = False
                     return overall_time+sim_time,is_done

             if detection_pop>0 and sum(self.get_pops())>=detection_pop:
                 is_done = True
                 return overall_time+tot_t,is_done

        
         return overall_time+tot_t,is_done
