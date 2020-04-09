import config
import fileLoader as FL
import operator
import matplotlib.pyplot as plt


def getAllTargetClassesSorted(allTargets):
    notInFileCounter = 0
    allTargetClasses = FL.loadAndBuildDictionary(config.TARGET_CLASSES_PATH, 4, "|")
    targetClasses = {}
    for target in allTargets:
        if target in allTargetClasses:
            splitClasses = allTargetClasses[target].split(".")
            domainClass = splitClasses[0]
            if domainClass in targetClasses:
                targetClasses[domainClass] += 1
            else:
                if domainClass == "Protein Kinase Superfamily":
                    print(target)
                targetClasses[domainClass] = 1
        else:
            notInFileCounter += 1

    print("Targets not found in file: " + str(notInFileCounter))
    return sorted(targetClasses.items(), key=operator.itemgetter(1), reverse=True)


def writeClassesToFile(path, targetClasses):
    file = open(path, "w")
    file.write("target class|generalized target class\n")
    for tClass, amount in targetClasses:
        if amount >= 10 and tClass != '':
            file.write(tClass + "|" + "UNKNOWN\n")
    file.close()


def printTargetClasses(targetClasses):
    for entry in targetClasses:
        print(entry)

    limitedAmount = 200
    counter = 0
    otherValues = 0
    limitedTargetClasses = {}
    for key, value in targetClasses:
        if key != '' and 3 <= value:
            limitedTargetClasses[key] = value
            counter += 1
        else:
            otherValues += value

    limitedTargetClasses['Other'] = otherValues
    print(limitedTargetClasses)
    plt.pie([v for v in limitedTargetClasses.values()], labels=[k for k in limitedTargetClasses.keys()], autopct=None)
    plt.show()


allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
targetClassesSorted = getAllTargetClassesSorted(allTargets)
#writeClassesToFile(config.TARGET_CLASSES_GENERAL_GROUPS, targetClassesSorted)
printTargetClasses(targetClassesSorted)
