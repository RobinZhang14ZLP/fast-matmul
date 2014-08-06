# Generate all of the algorithms needed for the comparison of addition performance.
# Usage:
#
#           ./gen_adds_perf_algs.sh


# Where the U, V, W versions of the algorithms live
ALGS_DIR=algorithms
ELIM_DIR=${ALGS_DIR}/eliminated

# Where we are going to place the generated code.
OUT_DIR=../algorithms/special
mkdir -p ${OUT_DIR}

SCRIPT=gen.py
ALG1=fast424_26_257
ALG2=fast433_29_234
# No common subexpression elimination.
# Write-once
OPTIONS=0
python ${SCRIPT} ${ALGS_DIR}/grey424-26-257 4,2,4 ${OUT_DIR}/${ALG1}_wo.hpp ${OPTIONS} ${ALG1}_wo
python ${SCRIPT} ${ALGS_DIR}/grey433-29-234 4,3,3 ${OUT_DIR}/${ALG2}_wo.hpp ${OPTIONS} ${ALG2}_wo
# Streaming
OPTIONS=1
python ${SCRIPT} ${ALGS_DIR}/grey424-26-257 4,2,4 ${OUT_DIR}/${ALG1}_st.hpp ${OPTIONS} ${ALG1}_st
python ${SCRIPT} ${ALGS_DIR}/grey433-29-234 4,3,3 ${OUT_DIR}/${ALG2}_st.hpp ${OPTIONS} ${ALG2}_st
# Pairwise
OPTIONS=2
python ${SCRIPT} ${ALGS_DIR}/grey424-26-257 4,2,4 ${OUT_DIR}/${ALG1}_pw.hpp ${OPTIONS} ${ALG1}_pw
python ${SCRIPT} ${ALGS_DIR}/grey433-29-234 4,3,3 ${OUT_DIR}/${ALG2}_pw.hpp ${OPTIONS} ${ALG2}_pw

ALG1=fast424_26_206
ALG2=fast433_29_195
# Yes common subexpression elimination.
# Write-once
OPTIONS=0
python ${SCRIPT} ${ELIM_DIR}/grey424-26-206 4,2,4 ${OUT_DIR}/${ALG1}_wo.hpp ${OPTIONS} ${ALG1}_wo
python ${SCRIPT} ${ELIM_DIR}/grey433-29-195 4,3,3 ${OUT_DIR}/${ALG2}_wo.hpp ${OPTIONS} ${ALG2}_wo
# Streaming
OPTIONS=1
python ${SCRIPT} ${ELIM_DIR}/grey424-26-206 4,2,4 ${OUT_DIR}/${ALG1}_st.hpp ${OPTIONS} ${ALG1}_st
python ${SCRIPT} ${ELIM_DIR}/grey433-29-195 4,3,3 ${OUT_DIR}/${ALG2}_st.hpp ${OPTIONS} ${ALG2}_st
# Pairwise
OPTIONS=2
python ${SCRIPT} ${ELIM_DIR}/grey424-26-206 4,2,4 ${OUT_DIR}/${ALG1}_pw.hpp ${OPTIONS} ${ALG1}_pw
python ${SCRIPT} ${ELIM_DIR}/grey433-29-195 4,3,3 ${OUT_DIR}/${ALG2}_pw.hpp ${OPTIONS} ${ALG2}_pw

