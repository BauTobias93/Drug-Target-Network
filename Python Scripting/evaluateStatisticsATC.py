
import config
import fileLoader as FL
import operator
import matplotlib.pyplot as plt

def printAllATC(drugs):
    categories = {}
    counterNoATC = 0
    for d in drugs:
        if not drugs[d].atc:
            counterNoATC += 1
        else:
            for atc in drugs[d].atc:
                if atc[0:1] not in categories:
                    categories[atc[0:1]] = 1
                else:
                    categories[atc[0:1]] += 1

    #categories = sorted(categories.items(), key=operator.itemgetter(1),reverse=True)
    print(categories)
    plt.pie([v for v in categories.values()], labels=[k for k in categories.keys()], autopct=None)
    plt.show()




drugs = FL.getAllDrugs(config.DRUG_FILE_PATH)
allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
printAllATC(drugs)