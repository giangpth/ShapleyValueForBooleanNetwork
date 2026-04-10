import timeit
from . import BooleanShapleyTest, BooleanShapleyAnalysis, testEachModel, runBatchTest

if __name__ == "__main__":
    start = timeit.default_timer()
    # BooleanShapleyAnalysis()
    BooleanShapleyTest()
    # runBatchTest('realbenchmark.txt', 'cycle1_1_nolog.json')
    # testEachModel('data/realbenchmark/1.txt', 'Apoptosis')
    print("-----------TOTAL RUNNING TIME---------------")
    print(timeit.default_timer() - start)  