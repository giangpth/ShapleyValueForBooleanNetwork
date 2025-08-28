import timeit
from . import BooleanShapleyAnalysis 

if __name__ == "__main__":
    start = timeit.default_timer()
    BooleanShapleyAnalysis()
    print("-----------RUNNING TIME---------------")
    print(timeit.default_timer() - start)