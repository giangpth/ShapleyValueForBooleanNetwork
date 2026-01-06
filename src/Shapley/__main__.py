import timeit
from . import BooleanShapleyTest, BooleanShapleyAnalysis, testEachModel, percentageTest

if __name__ == "__main__":
    start = timeit.default_timer()
    # BooleanShapleyAnalysis()
    # BooleanShapleyTest()
    percentageTest('realbenchmark_2.txt', 'output_2.json')
    # testEachModel('data/realbenchmark/1.txt', 'Apoptosis')
    print("-----------TOTAL RUNNING TIME---------------")
    print(timeit.default_timer() - start)  