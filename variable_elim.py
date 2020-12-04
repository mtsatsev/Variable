"""
@Author: Joris van Vugt, Moira Berens, Leonieke van den Bulk

Class for the implementation of the variable elimination algorithm.
http://www.nbertagnolli.com/jekyll/update/2016/05/23/Bayes_Nets.html
https://github.com/RaptorMai/bayesian-network-variable-elimination-gibbs-sampling/blob/master/Factor.py

"""

from Factor import product,marginalize,reduce
from read_bayesnet import BayesNet
import sys

class VariableElimination():

    def __init__(self, network):
        """
        Initialize the variable elimination algorithm with the specified network.
        Add more initializations if necessary.

        """
        self.network = network


    def elimination_order(self, query, function):
        if query not in self.network.nodes:
            raise ValueError("The query variable {} is not in the network".format(query))

        order = [x for x in network.nodes if x != query]
        if not function:
            return order
        else:
            order.sort(key=function,reverse=True)
            return order

    def fewest_factors(self):
        def function(node):
            result = 0
            for n in network.nodes:
                factors = network.probabilities[n].columns.tolist()[:-1]
                if node in factors:
                    result+=1
            return result
        return function

    def least_arcs(self):
        def function(node):
            result = network.parents[node]
            return result
        return function

    def run(self, query, observed, elim_order):
        """
        Use the variable elimination algorithm to find out the probability
        distribution of the query variable given the observed variables

        Input:
            query:      The query variable
            observed:   A dictionary of the observed variables {variable: value}
            elim_order: Either a list specifying the elimination ordering
                        or a function that will determine an elimination ordering
                        given the network during the run

        Output: A variable holding the probability distribution
                for the query variable

        """
        out = sys.stdout
        logs = open("notes.txt", "w")
        sys.stdout = logs
        print("++++++++++++++++++++++ The nodes we have ++++++++++++++++++++++")
        print("++++++++++++++++++++++"  + network.nodes.__str__() + "++++++++++++++++++++++")
        print("++++++++++++++++++++++ The query variable is: "+query+" ++++++++++++++++++++++")

        order = [x for x in elim_order if x != query]
        print("++++++++++++++++++++++ Our elimination order is: " + str(elim_order)+"\n")
        facs = []
        print("++++++++++++++++++++++ Remove the observed variable(s) which is/are: "+str(observed))
        for node in observed:
            self.observe(self.network, node, observed[node])
            self.removeFromNet(self.network, node)
            if node in order:
                order.remove(node)
                xs = [x for x in self.network.nodes if node != x]
                facs = facs + (xs)

        facs = set(facs)
        formulas = [x for x in facs]
        elim_order = order
        if len(formulas) == 0:
            formulas =[x for x in self.network.nodes]
        
        logs.write("++++++++++++++++++++++We are going to be working with++++++++++++++++++++++\n ")

        self.properPrint(self.network.probabilities)
        while (elim_order):
            node = elim_order.pop(0)
            logs.write("++++++++++++++++++++++ We are eliminating" + node + " ++++++++++++++++++++++\n")
            factor, formulas = self.multiply(node, formulas)
            logs.write("++++++++++++++++++++++ Now we are left with ++++++++++++++++++++++\n")

            logs.write("++++++++++++++++++++++" + str(formulas) + " ++++++++++++++++++++++\n")
            logs.write("++++++++++++++++++++++ Now We are working with ++++++++++++++++++++++\n ")
            self.properWrite(formulas)

            if factor != None:
                logs.write("++++++++++++++++++++++ Marginalize the table ++++++++++++++++++++++\n")
                self.network.probabilities[factor] = marginalize(node, self.network.probabilities[factor])
                formulas = formulas + [factor]
            logs.write("++++++++++++++++++++++Now We are working with++++++++++++++++++++++\n ")
            self.properWrite(formulas)

        for fac in formulas:
            if self.network.probabilities[fac].shape[0] == 1:
                formulas.remove(fac)
        logs.write("++++++++++++++++++++++Before the end++++++++++++++++++++++\n ")
        self.properWrite(formulas)
        result = formulas[0]
        for r in range(1, len(formulas)):
            self.network.probabilities[result] = product(self.network.probabilities[result], self.network.probabilities[formulas[r]])

        sum = 0
        for prob in network.probabilities[result]['prob']:
            sum += prob
        network.probabilities[result]['prob'] = network.probabilities[result]['prob']/sum
        logs.write("++++++++++++++++++++++ The result is: ++++++++++++++++++++++\n ")
        print(network.probabilities[result])
        sys.stdout = out
        logs.close()

        return result

    def observe(self,network,node,evidence):
        for n in network.nodes:
            if node in network.probabilities[n].columns.tolist():
                network.probabilities[n] = reduce(node,network.probabilities[n],evidence)
        network.probabilities.pop(node)

    def multiply(self,node, formulas):
        first = True
        newFactor = None
        for factor in formulas:
            if node in network.probabilities[factor].columns.tolist():
                if first:
                    first = False
                    newFactor = factor
                else:
                    network.probabilities[newFactor] = product(network.probabilities[newFactor],
                                                               network.probabilities[factor])

                formulas[formulas.index(factor)] = None

        while None in formulas:
            formulas.remove(None)
        return newFactor, formulas

    def removeFromNet(self,network,node):
        for probs in network.probabilities:
            if node in network.probabilities[probs]:
                network.probabilities[probs] = network.probabilities[probs].drop(node,axis=1)

    def properPrint(self,probabilities):
        for prob in probabilities:
            print(prob)
            print(network.probabilities[prob])

    def properWrite(self,formulas):
       for factor in formulas:
           print(self.network.probabilities[factor])

network = BayesNet('cancer.bif')
VE = VariableElimination(network)
query = "Cancer"
elim_order = VE.elimination_order(query,VE.fewest_factors())
observed = {"Pollution":"low"}
result = VE.run(query,observed,elim_order)
print("Result")
print(VE.network.probabilities[result])


