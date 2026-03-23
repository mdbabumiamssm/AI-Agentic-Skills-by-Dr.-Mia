# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# get the plays and measure the performance of query

import pymysql
import time

myConnection = pymysql.connect(
    host="localhost", user="root", password="root", db="shakespeare")

cur = myConnection.cursor()
start_time = time.time()

# Part 1: Get the plays from database
cur.execute("SELECT play_text FROM amnd;")

for line in cur.fetchall():
    print(line[0])

end_time = time.time()

# Part 2: process for Query Performance Calculation
cur.execute('SELECT COUNT(line_number) FROM amnd;')
numPlayLines = cur.fetchall()[0][0]
print(numPlayLines, 'rows')

# calculate query execution time
queryExecTime = end_time - start_time
print("Total query time: ", queryExecTime)

queryTimePerLine = queryExecTime / numPlayLines
print("Query time per line: ", queryTimePerLine)

# record query execution time into performance table
insertPerformanceSQL = "INSERT INTO performance VALUES('READ',%s);"
cur.execute(insertPerformanceSQL, queryTimePerLine)

myConnection.commit()
myConnection.close()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
