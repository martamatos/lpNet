library(RUnit)
source("../scripts/build_timeSeriesData.r")
source("../src/sourceDir.R")
sourceDir("../src")

# directory where test files are located
dirs = c("/home/dx/Thesis/lpNet_Code/tests")
testFileRegexp = "test_bootstrap.r"
testFuncRegexp = "test*"


testSuite = defineTestSuite("test lpNet", dirs, testFileRegexp, testFuncRegexp)
testResult = runTestSuite(testSuite, useOwnErrorHandler = TRUE)

printTextProtocol(testResult, showDetails = TRUE)

#verbose = getOption("RUnit")$verbose)



