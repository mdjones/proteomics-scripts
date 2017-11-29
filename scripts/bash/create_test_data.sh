#!/bin/bash


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TEST_DATA_DIR=${DIR}/../../test/data/
PYTHON=/usr/bin/python
ANALYZE_QUANT_SCRIPT="${DIR}/../../data/NomuraProtScript11162017/analyze_quantCompare.py"
PEP_LIST_IN="${DIR}/../../data/NomuraProtScript11162017/peptideList_EN80.csv"


## Create test data
grep -E -f ${DIR}/peptideListGrepPatterns.txt ${PEP_LIST_IN} > ${TEST_DATA_DIR}/peptideList.csv

CMD="${PYTHON} ${ANALYZE_QUANT_SCRIPT} ${TEST_DATA_DIR}/peptideList.csv ${TEST_DATA_DIR}/results.csv"
echo ${CMD}
${CMD}

CMD="${PYTHON} ${ANALYZE_QUANT_SCRIPT} -v ${TEST_DATA_DIR}/peptideList.csv ${TEST_DATA_DIR}/resultsverbose.csv"
echo ${CMD}
${CMD}

