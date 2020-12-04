import pandas as pd
import numpy as np

def marginalize(node,probs):
    print(probs.shape)
    print(node)
    print("HEYLYA")
    print(probs)
    if probs.shape == (2,2):
        return probs
    probability = probs['prob']
    cpt=probs.drop([node,'prob'],axis=1)
    marginal = pd.DataFrame(columns=cpt.columns.tolist())
    probs = []
    length = cpt.shape[1]
    while cpt.shape[0] > 0:
        positions = [x for x in range(0,cpt.shape[0]) if sum(cpt.iloc[0]==cpt.iloc[x]) == cpt.shape[1]]

        probs.append(sum(probability[probability.index[positions]]))
        marginal = marginal.append(cpt[:1])

        cpt=cpt.drop(cpt.index[positions],axis=0)
        probability=probability.drop(probability.index[positions],axis=0)

    marginal.insert(length,'prob',probs)
    return marginal


def product(probs1, probs2):
    intersec = list(filter(lambda x: x != 'prob', np.intersect1d(probs1.columns, probs2.columns)))

    probone = probs1['prob']
    probtwo = probs2['prob']

    cpt1 = probs1.drop('prob', axis=1)
    cpt2 = probs2.drop('prob', axis=1)

    prod = cpt2.join(cpt1.set_index(intersec), on=intersec)
    probs = []

    for i in range(0, cpt1.shape[0]):
        for j in range(0, cpt2.shape[0]):
            if sum(cpt1[intersec].iloc[i] == cpt2[intersec].iloc[j]) == len(intersec):
                probs.append(probone.iloc[i] * probtwo.iloc[j])

    prod['prob'] = probs
    return prod

def reduce(node,probs,evidence=None):
    if not evidence:
        return probs
    columns = probs.columns.tolist()
    try:
        i = columns.index(node)
        reduced = pd.DataFrame(list(filter(lambda x: x[i] == evidence,probs.values)),columns=columns)
        return reduced
    except:
        print("{} is not found in the probabilities".format(node))