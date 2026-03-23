# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Delete Stage Directions in Plays text such as "Exit", "Enter", etc

import pymysql
import time

myConnection = pymysql.connect(
    host="localhost", user="root", password="root", db="shakespeare")

cur = myConnection.cursor()
start_time = time.time()

cur.execute('SELECT COUNT(line_number) FROM amnd;')
numPlayLines_Before_Delete = cur.fetchall()[0][0]


# delete stage directions
# RLIKE : regular expression like, ^ means starts with , | means OR
cur.execute(
    "DELETE FROM amnd WHERE play_text RLIKE '^enter|^exit|^act|^scence|^exeunt';")

print("Deleting lines..")

end_time = time.time()

myConnection.commit()


cur.execute('SELECT COUNT(line_number) FROM amnd;')
numPlayLines_After_Delete = cur.fetchall()[0][0]
numPlayLInes_Deleted = numPlayLines_Before_Delete - numPlayLines_After_Delete
print(numPlayLInes_Deleted, 'rows')

# calculate query execution time
queryExecTime = end_time - start_time
print("Total query time: ", queryExecTime)

queryTimePerLine = queryExecTime / numPlayLInes_Deleted
print("Query time per line: ", queryTimePerLine)

# record query execution time into performance table
insertPerformanceSQL = "INSERT INTO performance VALUES('DELETE',%s);"
cur.execute(insertPerformanceSQL, queryTimePerLine)

myConnection.commit()
myConnection.close()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
