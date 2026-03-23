# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# UK Road Safty 2015 data
import pymysql

myConnection = pymysql.connect(
    host="localhost", user="root", password="root", db="accidents")

cur = myConnection.cursor()

cur.execute(
    "SELECT vehicle_type FROM vehicle_types WHERE vehicle_type LIKE '%torcycle%';")

cycle_list = cur.fetchall()


selectSQL = ('''
SELECT vt.vehicle_type, a.accident_severity
FROM accident a
JOIN vehicles v ON a.accident_index = v.accident_index
JOIN vehicle_types vt ON v.vehicle_type = vt.vehicle_code
WHERE vt.vehicle_type LIKE %s
ORDER BY a.accident_severity;
''')


insert_SQL = ('''INSERT INTO accidents_median
VALUES(%s, %s);''')


for cycle in cycle_list:
    cur.execute(selectSQL, cycle[0])
    accidents = cur.fetchall()

    # calculate median severity
    # divide the length of accidents /2 to find the median of accdients list
    quotient, remainder = divmod(len(accidents), 2)

    if remainder:
        # meaning odds number of items in accidents list
        median_severity = accidents[quotient][1]
    else:
        # even numbers of items in accidents list
        median_severity = (accidents[quotient]
                           [1] + accidents[quotient + 2][1]) / 2

    print("finding Median Severity for ", cycle[0])

    # insert the calculated median severity into table
    cur.execute(insert_SQL, (cycle[0], median_severity))

myConnection.commit()
myConnection.close()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
