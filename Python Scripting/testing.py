
import numpy as np
import config



matrixTargetGODistance = np.load(config.MATRIX_TARGET_DISTANCE)
np.savetxt("./test/" + config.MATRIX_TARGET_DISTANCE[2:-4] + ".csv", matrixTargetGODistance, delimiter=";", fmt=('%.3f'))




